#!/usr/bin/env bash
set -euo pipefail

AUTHOR="Codruta Toadere"
EMAIL="codrutatoadere@gmail.com"

init()
{
    use_locale=$1
    echo "..: Creating ${use_locale} translation."

    lang=${use_locale:0:2}
    msginit --locale=${use_locale} --input=messages.pot --output=messages-${lang}.po
}

compile()
{
    use_locale=$1
    echo "..: Building ${use_locale} translation binaries"

    lang=${use_locale:0:2}
    msgmerge messages.pot \
        --lang=${use_locale} \
        --sort-output \
        --verbose \
        --update "messages-${lang}.po"
    mkdir -p "./${lang}/LC_MESSAGES"
    msgfmt "messages-${lang}.po" \
        --verbose \
        --output="./${lang}/LC_MESSAGES/fbg-simulation-pyqt.mo"
}

# collect all translatable strings
xgettext gui/*.py \
    --sort-output \
    --foreign-user \
    --copyright-holder="${AUTHOR}" \
    --package-name="fbg-simulation-pyqt" \
    --package-version="v2.0" \
    --msgid-bugs-address="${EMAIL}" \
    --verbose \
    --output="./translation/messages.pot"

# cd ./translation
# for language in ro_RO en_US; do
#     compile $language
# done
# cd ..
