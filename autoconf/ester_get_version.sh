#!/bin/sh

versionfile=$1

do_subst=` sed \
    -e 's/^major/ESTER_MAJOR_VERSION/' \
    -e 's/^minor/ESTER_MINOR_VERSION/' \
    -e 's/^release/ESTER_RELEASE_VERSION/' \
    -e 's/^greek/ESTER_GREEK_VERSION/' < $versionfile`

eval "$do_subst"

ESTER_VERSION="${ESTER_MAJOR_VERSION}.${ESTER_MINOR_VERSION}.${ESTER_RELEASE_VERSION}${ESTER_GREEK_VERSION}"

echo $ESTER_VERSION
