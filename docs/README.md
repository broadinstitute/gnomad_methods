# Documentation

Documentation for gnomad_hail is generated using [Sphinx](https://www.sphinx-doc.org/en/master/).

To build the documentation, run:

```
cd /path/to/gnomad_hail
pip install -r requirements.txt
pip install -r docs/requirements.txt
./docs/build.sh
```

The generated HTML is placed in docs/html. To view the documentation, run:

```
python3 -m http.server -d docs/html
```

## References

- [Sphinx documentation](https://www.sphinx-doc.org/en/master/)
- [Sphinx reStructuredText documentation](https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html)
- [docutils reStructuredText documentation](https://docutils.sourceforge.io/rst.html)
