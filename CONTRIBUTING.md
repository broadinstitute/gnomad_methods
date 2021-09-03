# Contributing

## Setting up a development environment

- Install development dependencies.

  ```
  python3 -m pip install -r requirements-dev.txt
  ```

- Install [pre-commit](https://pre-commit.com/) hooks.

  ```
  python3 -m pre_commit install
  ```

## Running Pylint

Use [Pylint](https://www.pylint.org/) to check for some types of errors.

```
./lint
```

To disable warnings, use:

```
./lint --disable=W
```

## Running Black

Use [Black](https://black.readthedocs.io/) to format code.

```
black gnomad
```

## Running pydocstyle

Use [pydocstyle](https://www.pydocstyle.org/) to check that docstrings conform to [PEP 257](https://www.python.org/dev/peps/pep-0257/).

```
pydocstyle gnomad
```

## Running tests

Run tests using [pytest](https://docs.pytest.org/en/stable/).

```
python -m pytest
```

## Building documentation

See instructions in [docs/README.md](./docs/README.md).

## Publishing a release

https://packaging.python.org/guides/distributing-packages-using-setuptools/#packaging-your-project

- Make sure that the [changelog](./CHANGELOG.md) is up to date.

  To see changes since the last release, use:

  ```
  LAST_RELEASE_TAG=$(git tag --list --sort=-committerdate | head -n1)
  git log $LAST_RELEASE_TAG..
  ```

- Update version in setup.py, replace "unreleased" heading in changelog with the version number, and commit.
  Push changes to the master branch on GitHub or submit a pull request.
  The new version number should be based on changes since the last release.

  https://semver.org/

  https://packaging.python.org/guides/distributing-packages-using-setuptools/#semantic-versioning-preferred

- Once the version has been updated in the master branch on GitHub, tag the release.
  The version tag should be applied to a commit on the master branch.

  The version number in the tag must match the version number in setup.py.

  ```
  git checkout master
  git pull

  git tag v<version>
  git push origin v<version>
  ```

  When a `v<VERSION_NUMBER>` tag is pushed to GitHub, an Actions workflow will automatically publish the tagged code to PyPI.

  New releases will automatically be posted in the #gnomad_notifications Slack channel (via the RSS Slack app).

### Manually publishing a release

- Install tools.

  ```
  python -m pip install --upgrade setuptools wheel twine
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
