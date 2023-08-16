#!/usr/bin/env bash
set -euo pipefail

AUTHOR="Codruta Toadere"
EMAIL="codrutatoadere@gmail.com"

init()
{
    msginit --locale=en_US --input=messages.pot --output=messages-en.po
    msginit --locale=ro_RO --input=messages.pot --output=messages-ro.po
}

compile()
{
    language_code=$1
    lang=${language_code:0:2}
    msgmerge messages.pot\
        --lang=${language_code} \
        --sort-output \
        --update "messages-${lang}.po"

    mkdir -p "./${lang}/LC_MESSAGES"
    msgfmt "messages-${lang}.po" --output="./${lang}/LC_MESSAGES/fbg-simulation-pyqt.mo"
}

cd ./translation

# collect all translatable strings
xgettext ../gui/*.py \
    --sort-output \
    --foreign-user \
    --copyright-holder="${AUTHOR}" \
    --package-name="fbg-simulation-pyqt" \
    --package-version="v2.0" \
    --msgid-bugs-address="${EMAIL}" \
    --output="messages.pot" \

for language in ro_RO en_US; do
    compile $language
done

cd ..
