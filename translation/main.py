import gettext


_translators = dict()


def install(lang: str):
    gettext.bindtextdomain("fbg-simulation-pyqt", "./translation")
    gettext.textdomain("fbg-simulation-pyqt")

    global _translators
    _translators = dict(
        ro=gettext.translation("fbg-simulation-pyqt", "./translation", languages=["ro"]),
        en=gettext.translation("fbg-simulation-pyqt", "./translation", languages=["en"]),
    )
    _translators.get(lang).install()
