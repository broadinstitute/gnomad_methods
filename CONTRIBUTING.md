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

## Running isort

Use [isort](https://pycqa.github.io/isort/) to format imports.

```
isort gnomad
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

- Update version in setup.py and commit.
  Push changes to the main branch on GitHub or submit a pull request.
  The new version number should be based on changes since the last release.

  https://semver.org/

  https://packaging.python.org/guides/distributing-packages-using-setuptools/#semantic-versioning-preferred

- Once the version has been updated in the main branch on GitHub, tag the release.
  The version tag should be applied to a commit on the main branch.

  The version number in the tag must match the version number in setup.py.

  The tag can be created using GitHub by [creating a GitHub Release with a new tag](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository#creating-a-release).

  Alternatively, the tag can be created locally and pushed to GitHub using:

  ```
  git checkout main
  git pull

  git tag v<version>
  git push origin v<version>
  ```

- When a `v<VERSION_NUMBER>` tag is pushed to GitHub, an Actions workflow will automatically publish the tagged code to PyPI.

  New releases will automatically be posted in the #gnomad_notifications Slack channel (via the RSS Slack app).

- If the version tag was created locally and pushed to GitHub, then a [GitHub Release](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases) will be automatically created with
  [release notes generated from pull requests](https://docs.github.com/en/repositories/releasing-projects-on-github/automatically-generated-release-notes).

  Check the [Releases page](https://github.com/broadinstitute/gnomad_methods/releases) to make sure the generated release notes look ok
  and edit them if necessary.

  If needed, to see commits since the last release, use:

  ```
  LAST_RELEASE_TAG=$(git tag --list --sort=-committerdate | head -n1)
  git log $LAST_RELEASE_TAG..
  ```

  Especially if there are many changes in the release, organize the changelog by type.

  See https://keepachangelog.com/en/1.0.0/#how for more information.

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
