import io


def get_slack_info():
    import getpass
    from slackclient import SlackClient
    # import os
    try:
        from gnomad.slack_creds import slack_token  # pylint: disable=import-error,no-name-in-module
    except Exception:
        return None

    # slack_token = os.environ["SLACK_API_TOKEN"]
    sc = SlackClient(slack_token)
    user = getpass.getuser()
    if user.startswith('konrad'): user = 'konradjk'
    users = [x['name'] for x in sc.api_call("users.list")['members']]
    default_channel = '#gnomad' if user not in users else '@' + user
    return sc, default_channel


def get_slack_channel_id(sc, channel):
    channel_id = [x for x in sc.api_call('channels.list')['channels'] if x['name'] == channel]
    return channel_id[0]['id'] if len(channel_id) and 'id' in channel_id[0] else None


def get_slack_user_id(sc, user):
    channel_id = [x for x in sc.api_call('users.list')['members'] if x['name'] == user]
    return channel_id[0]['id'] if len(channel_id) and 'id' in channel_id[0] else None


def send_message(channels=None, message="Your job is done!", icon_emoji=':woohoo:'):
    sc, default_channel = get_slack_info()

    if not isinstance(channels, list):
        channels = [channels]
    for channel in channels:
        sc.api_call(
            "chat.postMessage",
            channel=channel,
            text=message,
            icon_emoji=icon_emoji,
            parse='full'
        )

        
def send_snippet(channels=None, content='', filename='data.txt'):
    sc, default_channel = get_slack_info()

    if isinstance(channels, str):
        channels = [channels]
    elif channels is None:
        channels = [default_channel]

    for channel in channels:
        if channel.startswith('@'):
            channel_id = get_slack_user_id(sc, channel.lstrip('@'))
        else:
            channel_id = get_slack_channel_id(sc, channel.lstrip('#'))

        try:
            return sc.api_call("files.upload",
                               channels=channel_id,
                               content=content,
                               filename=filename)
        except Exception:
            print('Slack connection fail. Was going to send:')
            print(content)


def send_file(p, channel, filename='plot.png', tfile_name='/home/hail/plot.png'):
    from bokeh.io import export_png
    sc, default_channel = get_slack_info()
    if channel.startswith('@'):
        channel_id = get_slack_user_id(sc, channel.lstrip('@'))
    else:
        channel_id = get_slack_channel_id(sc, channel.lstrip('#'))

    export_png(p, filename=tfile_name)
    with open(tfile_name, 'rb') as f:
        sc.api_call("files.upload",
                    channels=channel_id,
                    file=io.BytesIO(f.read()),
                    filename=filename
                    )


def try_slack(target, func, *args):
    import sys
    import os
    import traceback
    import time
    process = os.path.basename(sys.argv[0])
    try:
        func(*args)
        send_message(target, 'Success! {} finished!'.format(process))
    except Exception as e:
        try:
            emoji = ':white_frowning_face:'
            if len(str(e)) > 3500:  # Slack message length limit (from https://api.slack.com/methods/chat.postMessage)
                filename = 'error_{}_{}.txt'.format(process, time.strftime("%Y-%m-%d_%H:%M"))
                snippet = send_snippet(target, traceback.format_exc(), filename=filename)
                if 'file' in snippet:
                    if 'SparkContext was shut down' in str(e) or 'connect to the Java server' in str(e):
                        send_message(target, 'Job ({}) cancelled - see {} for error log'.format(process, snippet['file']['url_private']), ':beaker:')
                    else:
                        send_message(target, 'Job ({}) failed :white_frowning_face: - see {} for error log'.format(process, snippet['file']['url_private']), emoji)
                else:
                    send_message(target, 'Snippet failed to upload: {}'.format(snippet), emoji)
            else:
                send_message(target, 'Job ({}) failed :white_frowning_face:\n```{}```'.format(process, traceback.format_exc()), emoji)
            raise e
        except ImportError as f:
            print("ERROR: missing slackclient. But here's the original error:")
            raise e