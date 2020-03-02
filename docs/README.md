# Documentation

Documentation for gnomad is generated using [Sphinx](https://www.sphinx-doc.org/en/master/).

To build the documentation, run:

```
cd /path/to/gnomad
pip install -r requirements.txt
pip install -r docs/requirements.txt
./docs/build.sh
```

The generated HTML is placed in docs/html. To view the documentation, run:

```
(cd docs/html && python3 -m http.server)
```

and open http://localhost:8000 in a web browser.

## References

- [Sphinx documentation](https://www.sphinx-doc.org/en/master/)
- [Sphinx reStructuredText documentation](https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html)
- [docutils reStructuredText documentation](https://docutils.sourceforge.io/rst.html)
