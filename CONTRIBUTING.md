# Contributing

## Running Pylint

Use [Pylint](https://www.pylint.org/) to check for some types of errors.

```
python -m pip install pylint
./lint
```

To disable warnings, use:

```
./lint --disable=W
```

## Building documentation

See instructions in [docs/README.md](./docs/README.md).

## Publishing a release

https://packaging.python.org/guides/distributing-packages-using-setuptools/#packaging-your-project

- Install tools.

  ```
  python -m pip install --upgrade setuptools wheel twine
  ```

- Update version in setup.py and commit.
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
