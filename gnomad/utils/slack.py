import os
import sys
import time
import traceback
import typing
from contextlib import contextmanager


from slack import WebClient


class SlackClient:
    """
    Slack API client.

    :param token: Slack API token
    """

    def __init__(self, token: str):
        self._client = WebClient(token=token)
        self._display_name_map = None

    def _load_display_name_map(self):
        display_name_map = {}
        response = self._client.users_list(limit=100)
        for user in response["members"]:
            if not (user["deleted"] or user["is_bot"]):
                display_name_map[user["profile"]["display_name"]] = user["id"]

        while response["response_metadata"]["next_cursor"]:
            next_cursor = response["response_metadata"]["next_cursor"]
            response = self._client.users_list(cursor=next_cursor, limit=100)
            for user in response["members"]:
                if not (user["deleted"] or user["is_bot"]):
                    display_name_map[user["profile"]["display_name"]] = user["id"]

        self._display_name_map = display_name_map

    def _get_direct_message_channel(self, user: str):
        if not self._display_name_map:
            self._load_display_name_map()

        if user.startswith("@"):
            user = user[1:]

        try:
            user_id = self._display_name_map[user]
        except KeyError:
            raise ValueError(f"User '{user}' not found in this workspace")
        else:
            response = self._client.conversations_open(users=[user_id])
            return response["channel"]["id"]

    def send_file(
        self,
        to: typing.Union[str, typing.Iterable[str]],
        file: typing.Optional[str] = None,
        content: typing.Optional[str] = None,
        filename: str = "data.txt",
        filetype: str = "text",
        comment: typing.Optional[str] = None,
    ):
        """
        Send a file to Slack channel(s) and/or user(s).

        :param to: Channel(s) (prefixed with '#') and/or user(s) (prefixed with '@') to send message to
        :param file: Path of file to upload
        :param content: File content to upload
        :param filename: Filename of file
        :param filetype: File type identifier
        :param comment: Text for message sharing file
        """
        if not (content or file) or (content and file):
            raise ValueError(
                "One, but not both, of 'content' or 'file' must be provided"
            )

        if isinstance(to, str):
            to = [to]

        for channel_or_user in to:
            if channel_or_user.startswith("@"):
                channel = self._get_direct_message_channel(channel_or_user)
            else:
                channel = channel_or_user

            optional_args = {}
            if file:
                optional_args["file"] = file
            else:
                optional_args["content"] = content

            if comment:
                optional_args["initial_comment"] = comment

            self._client.files_upload(
                channels=channel, filename=filename, filetype=filetype, **optional_args,
            )

    def send_message(
        self,
        to: typing.Union[str, typing.Iterable[str]],
        message: str,
        icon_emoji: typing.Optional[str] = None,
    ):
        """
        Send a message to Slack channel(s) and/or user(s).

        :param to: Channel(s) (prefixed with '#') and/or user(s) (prefixed with '@') to send message to
        :param message: Message content (long messages will be converted to snippets)
        :param icon_emoji: Emoji to use as icon for message
        """
        if isinstance(to, str):
            to = [to]

        for channel_or_user in to:
            if channel_or_user.startswith("@"):
                channel = self._get_direct_message_channel(channel_or_user)
            else:
                channel = channel_or_user

            if len(message) > 4000:
                self._client.files_upload(
                    channels=channel,
                    content=message,
                    filename="message.txt",
                    filetype="text",
                )
            else:
                optional_args = {}
                if icon_emoji:
                    optional_args["icon_emoji"] = icon_emoji

                self._client.chat_postMessage(
                    channel=channel, text=message, parse="full", **optional_args
                )


@contextmanager
def slack_notifications(token: str, to: typing.Union[str, typing.Iterable[str]]):
    """
    Send a Slack notification after some code runs.

    If the wrapped code block raises an exception, the notification will include the exception and stack trace.

    Example usage:

    .. code-block:: python

        with slack_notifications(token, "@username"):
            run_analysis()

    :param token: Slack API token
    :param to: Channel(s) (prefixed with '#') and/or user(s) (prefixed with '@') to send notification to
    """
    process = os.path.basename(sys.argv[0])
    try:
        yield

        slack_client = SlackClient(token)
        slack_client.send_message(
            to, f":white_check_mark: Success! {process} finished!"
        )
    except Exception as e:
        slack_client = SlackClient(token)
        slack_client.send_file(
            to,
            content=traceback.format_exc(),
            filename=f"error_{process}_{time.strftime('%Y-%m-%d_%H:%M')}.log",
            filetype="text",
            comment=f":x: Error in {process}",
        )

        raise
