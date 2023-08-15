from polib import pofile


_translations = None


def init_translations(lang_code="ro") -> None:
    po_db = pofile(f"./gui/messages-{lang_code}.po")

    global _translations
    if _translations is None:
        _translations = {entry.msgid: entry.msgstr for entry in po_db}
    else:
        raise RuntimeError("Translations engine already initialized")


def tr(msgid: str) -> str:
    if msgid not in _translations:
        print(f"WARNING: {msgid} is not found")
    return _translations.get(msgid, msgid)
