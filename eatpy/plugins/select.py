from typing import MutableMapping, Any, Union, Iterable
import fnmatch

from .. import shared


class Select(shared.Plugin):
    def __init__(
        self, include: Union[Iterable, str] = "*", exclude: Union[Iterable, str] = ()
    ):
        self.include = include if not isinstance(include, str) else (include,)
        self.exclude = exclude if not isinstance(exclude, str) else (exclude,)

    def initialize(self, variables: MutableMapping[str, Any], *args, **kwargs):
        eliminate = set()
        for name in variables:
            use = False
            for include_pattern in self.include:
                if fnmatch.fnmatch(name, include_pattern):
                    use = True
            for exclude_pattern in self.exclude:
                if fnmatch.fnmatch(name, exclude_pattern):
                    use = False
            if not use:
                eliminate.add(name)
        for name in eliminate:
            self.logger.info('Dropping %s from state presented to filter' % name)
            del variables[name]
