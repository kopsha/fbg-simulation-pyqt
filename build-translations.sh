#!/usr/bin/env bash
set -euxo pipefail

init()
{
    msginit --locale=en_US --input=messages.pot --output=messages-en.po
    msginit --locale=ro_RO --input=messages.pot --output=messages-ro.po
}

cd ./translation

# collect all translatable strings
xgettext --output=messages.pot ../gui/*.py

# update translation tables
msgmerge --lang=en_US --update messages-en.po messages.pot --sort-output
msgmerge --lang=ro_RO --update messages-ro.po messages.pot --sort-output

# compile translation binaries
mkdir -p ./en/LC_MESSAGES ./ro/LC_MESSAGES
msgfmt messages-en.po --output=./en/LC_MESSAGES/fbg-simulation-pyqt.mo
msgfmt messages-ro.po --output=./ro/LC_MESSAGES/fbg-simulation-pyqt.mo

cd ..
