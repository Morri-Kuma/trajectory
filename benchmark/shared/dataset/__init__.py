"""Dataset registries and preprocessors.

Importing this package registers all dataset and preprocessor subclasses.
"""

import importlib
import pkgutil

for _loader, _module_name, _is_pkg in pkgutil.walk_packages(
    __path__, prefix=f"{__name__}."
):
    importlib.import_module(_module_name)
