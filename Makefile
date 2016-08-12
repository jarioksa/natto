PKG = natto
R = R --slave -e
Rcmd = R CMD
## Version number
VERSION = $(shell R --slave -e 'cat(read.dcf("DESCRIPTION", fields="Version"))')

## for package, generate untracked dynamic files, build and clean

docs:
	$(R) 'devtools::document()'

build: docs
	cd ..;\
	$(Rcmd) build $(PKG)
	git clean -f

check: build
	cd ..;\
	$(Rcmd) check $(PKG)_$(VERSION).tar.gz

install: build
	cd ..;\
	$(Rcmd) INSTALL $(PKG)_$(VERSION).tar.gz
