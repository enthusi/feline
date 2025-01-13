# Minimal Makefile for Sphinx documentation

.PHONY: help clean html

# Directories
SOURCEDIR = .
BUILDDIR = _build

clean:
	@rm -rf $(BUILDDIR)

html:
	@sphinx-build -b html $(SOURCEDIR) $(BUILDDIR)/html
