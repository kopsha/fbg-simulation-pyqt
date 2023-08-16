import gettext


was_initialized = False
_ = gettext.gettext


def init_translations(lang_code="ro") -> None:
    gettext.bindtextdomain("fbg-simulation-pyqt", "./translation")
    gettext.textdomain("fbg-simulation-pyqt")
    lang = "en_US"
    gettext.translation("fbg-simulation-pyqt", "./translation", languages=[lang]).install()

    global was_initialized
    if was_initialized:
        raise RuntimeError("Translations engine already initialized")
    else:
        was_initialized = True
