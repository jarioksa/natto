## Version number
VERSION := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)

## for package, generate untracked dynamic files, build and clean

pkg:
	R -q -e 'devtools::document()'
	cd ..;\
	R CMD build natto
	git clean -f
