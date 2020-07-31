# Contributing

## Running Pylint

Use [Pylint](https://www.pylint.org/) to check for some types of errors.

```
python -m pip install -r requirements-dev.txt
./lint
```

To disable warnings, use:

```
./lint --disable=W
```

## Running Black

Use [Black](https://black.readthedocs.io/) to format code.

```
python -m pip install -r requirements-dev.txt
black gnomad
```

## Building documentation

See instructions in [docs/README.md](./docs/README.md).

## Publishing a release

https://packaging.python.org/guides/distributing-packages-using-setuptools/#packaging-your-project

- Install tools.

  ```
  python -m pip install --upgrade setuptools wheel twine
  ```

- Make sure that the [changelog](./CHANGELOG.md) is up to date.

  To see changes since the last release, use:

  ```
  LAST_RELEASE_TAG=$(git tag --list --sort=-committerdate | head -n1)
  git log $LAST_RELEASE_TAG..
  ```

- Update version in setup.py, replace "unreleased" heading in changelog with the version number, and commit.
  The new version number should be based on changes since the last release.

  https://semver.org/

  https://packaging.python.org/guides/distributing-packages-using-setuptools/#semantic-versioning-preferred

- Tag the release in git.

  ```
  git tag v<version>
  git push origin v<version>
  ```

- Remove any old builds.

  ```
  rm -r ./dist
  ```

- Build [source distribution](https://packaging.python.org/guides/distributing-packages-using-setuptools/#source-distributions)
  and [wheel](https://packaging.python.org/guides/distributing-packages-using-setuptools/#wheels).

  ```
  python setup.py sdist bdist_wheel
  ```

- Upload using [twine](https://pypi.org/project/twine/).

  ```
  twine upload dist/*
  ```

- Announce new release on Slack.
