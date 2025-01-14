# Minimal Makefile for Sphinx documentation

.PHONY: help clean html

# Directories
SOURCEDIR = .
BUILDDIR = _build


html:
	@rm -rf $(BUILDDIR)
	@sphinx-build -b html $(SOURCEDIR) $(BUILDDIR)/html
