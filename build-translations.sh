#!/usr/bin/env bash
set -euxo pipefail

# collect all translatable strings
xgettext --output=messages.pot ../gui/*.py
msginit --locale=en_US --input=messages.pot --output=messages-en.po
msginit --locale=ro_RO --input=messages.pot --output=messages-ro.po
