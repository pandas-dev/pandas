from functools import wraps
from copy import deepcopy
from traitlets import TraitError
from traitlets.config.loader import (
    Config,
)
from jupyter_core.application import JupyterApp
from jupyter_server.serverapp import ServerApp
from jupyter_server.extension.application import ExtensionApp
from .traits import NotebookAppTraits


def NBAPP_AND_SVAPP_SHIM_MSG(trait_name): return (
    "'{trait_name}' was found in both NotebookApp "
    "and ServerApp. This is likely a recent change. "
    "This config will only be set in NotebookApp. "
    "Please check if you should also config these traits in "
    "ServerApp for your purpose.".format(
        trait_name=trait_name,
    )
)


def NBAPP_TO_SVAPP_SHIM_MSG(trait_name): return (
    "'{trait_name}' has moved from NotebookApp to "
    "ServerApp. This config will be passed to ServerApp. "
    "Be sure to update your config before "
    "our next release.".format(
        trait_name=trait_name,
    )
)


def EXTAPP_AND_NBAPP_AND_SVAPP_SHIM_MSG(trait_name, extapp_name): return (
    "'{trait_name}' is found in {extapp_name}, NotebookApp, "
    "and ServerApp. This is a recent change. "
    "This config will only be set in {extapp_name}. "
    "Please check if you should also config these traits in "
    "NotebookApp and ServerApp for your purpose.".format(
        trait_name=trait_name,
        extapp_name=extapp_name
    )
)


def EXTAPP_AND_SVAPP_SHIM_MSG(trait_name, extapp_name): return (
    "'{trait_name}' is found in both {extapp_name} "
    "and ServerApp. This is a recent change. "
    "This config will only be set in {extapp_name}. "
    "Please check if you should also config these traits in "
    "ServerApp for your purpose.".format(
        trait_name=trait_name,
        extapp_name=extapp_name
    )
)


def EXTAPP_AND_NBAPP_SHIM_MSG(trait_name, extapp_name): return (
    "'{trait_name}' is found in both {extapp_name} "
    "and NotebookApp. This is a recent change. "
    "This config will only be set in {extapp_name}. "
    "Please check if you should also config these traits in "
    "NotebookApp for your purpose.".format(
        trait_name=trait_name,
        extapp_name=extapp_name
    )
)


def NOT_EXTAPP_NBAPP_AND_SVAPP_SHIM_MSG(trait_name, extapp_name): return (
    "'{trait_name}' is not found in {extapp_name}, but "
    "it was found in both NotebookApp "
    "and ServerApp. This is likely a recent change. "
    "This config will only be set in ServerApp. "
    "Please check if you should also config these traits in "
    "NotebookApp for your purpose.".format(
        trait_name=trait_name,
        extapp_name=extapp_name
    )
)


def EXTAPP_TO_SVAPP_SHIM_MSG(trait_name, extapp_name): return (
    "'{trait_name}' has moved from {extapp_name} to "
    "ServerApp. Be sure to update your config before "
    "our next release.".format(
        trait_name=trait_name,
        extapp_name=extapp_name
    )
)


def EXTAPP_TO_NBAPP_SHIM_MSG(trait_name, extapp_name): return (
    "'{trait_name}' has moved from {extapp_name} to "
    "NotebookApp. Be sure to update your config before "
    "our next release.".format(
        trait_name=trait_name,
        extapp_name=extapp_name
    )
)


# A tuple of traits that shouldn't be shimmed or throw any
# warnings of any kind.
IGNORED_TRAITS = ("open_browser", "log_level", "log_format", "default_url", "show_banner")


class NotebookConfigShimMixin:
    """A Mixin class for shimming configuration from
    NotebookApp to ServerApp. This class handles warnings, errors,
    etc.

    This class should be used during a transition period for apps
    that are switching from depending on NotebookApp to ServerApp.

    After one release cycle, this class can be safely removed
    from the inheriting class.

    TL;DR

    The entry point to shimming is at the `update_config` method.
    Once traits are loaded, before updating config across all
    configurable objects, this class injects a method to reroute
    traits to their *most logical* classes.

    This class raises warnings when:
        1. a trait has moved.
        2. a trait is redundant across classes.

    Redundant traits across multiple classes now must be
    configured separately, *or* removed from their old
    location to avoid this warning.

    For a longer description on how individual traits are handled,
    read the docstring under `shim_config_from_notebook_to_jupyter_server`.
    """

    @wraps(JupyterApp.update_config)
    def update_config(self, config):
        # Shim traits to handle transition from NotebookApp to ServerApp
        shimmed_config = self.shim_config_from_notebook_to_jupyter_server(
            config)
        super().update_config(shimmed_config)

    def shim_config_from_notebook_to_jupyter_server(self, config):
        """Reorganizes a config object to reroute traits to their expected destinations
        after the transition from NotebookApp to ServerApp.

        A detailed explanation of how traits are handled:

        1. If the argument is prefixed with `ServerApp`,
            pass this trait to `ServerApp`.
        2. If the argument is prefixed with `NotebookApp`,
            * If the argument is a trait of `NotebookApp` *and* `ServerApp`:
                1. Raise a warning—**for the extension developers**—that
                    there's redundant traits.
                2. Pass trait to `NotebookApp`.
            * If the argument is a trait of just `ServerApp` only
                (i.e. the trait moved from `NotebookApp` to `ServerApp`):
                1. Raise a "this trait has moved" **for the user**.
                3. Pass trait to `ServerApp`.
            * If the argument is a trait of `NotebookApp` only, pass trait
                to `NotebookApp`.
            * If the argument is not found in any object, raise a
                `"Trait not found."` error.
        3. If the argument is prefixed with `ExtensionApp`:
            * If the argument is a trait of `ExtensionApp`,
                `NotebookApp`, and `ServerApp`,
                1. Raise a warning about redundancy.
                2. Pass to the ExtensionApp
            * If the argument is a trait of `ExtensionApp` and `NotebookApp`,
                1. Raise a warning about redundancy.
                2. Pass to ExtensionApp.
            * If the argument is a trait of `ExtensionApp` and `ServerApp`,
                1. Raise a warning about redundancy.
                2. Pass to ExtensionApp.
            * If the argument is a trait of `ExtensionApp`.
                1. Pass to ExtensionApp.
            * If the argument is a trait of `NotebookApp` but not `ExtensionApp`,
                1. Raise a warning that trait has likely moved to NotebookApp.
                2. Pass to NotebookApp
            * If the arguent is a trait of `ServerApp` but not `ExtensionApp`,
                1. Raise a warning that the trait has likely moved to ServerApp.
                2. Pass to ServerApp.
            * else
                * Raise a TraitError: "trait not found."
        """
        extapp_name = self.__class__.__name__

        # Pop out the various configurable objects that we need to evaluate.
        nbapp_config = config.pop('NotebookApp', {})
        svapp_config = config.pop('ServerApp', {})
        extapp_config = config.pop(extapp_name, {})

        # Created shimmed configs.
        # Leave the rest of the config alone.
        config_shim = deepcopy(config)
        svapp_config_shim = {}
        nbapp_config_shim = {}
        extapp_config_shim = {}

        extapp_traits = (
            self.__class__.class_trait_names() +
            ExtensionApp.class_trait_names()
        )
        svapp_traits = ServerApp.class_trait_names()
        nbapp_traits = (
            NotebookAppTraits.class_trait_names() +
            ExtensionApp.class_trait_names()
        )

        # 1. Handle ServerApp traits.
        svapp_config_shim.update(svapp_config)

        # 2. Handle NotebookApp traits.
        warning_msg = None
        for trait_name, trait_value in nbapp_config.items():
            in_svapp = trait_name in svapp_traits
            in_nbapp = trait_name in nbapp_traits
            if trait_name in IGNORED_TRAITS:
                # Pass trait through without any warning message.
                nbapp_config_shim.update({trait_name: trait_value})
            elif in_svapp and in_nbapp:
                warning_msg = NBAPP_AND_SVAPP_SHIM_MSG(trait_name)
                nbapp_config_shim.update({trait_name: trait_value})
            elif in_svapp:
                warning_msg = NBAPP_TO_SVAPP_SHIM_MSG(trait_name)
                svapp_config_shim.update({trait_name: trait_value})
            elif in_nbapp:
                nbapp_config_shim.update({trait_name: trait_value})
            else:
                raise TraitError("Trait, {}, not found.".format(trait_name))

            # Raise a warning if it's given.
            if warning_msg:
                self.log.warning(warning_msg)

        # 3. Handle ExtensionApp traits.
        warning_msg = None
        for trait_name, trait_value in extapp_config.items():
            in_extapp = trait_name in extapp_traits
            in_svapp = trait_name in svapp_traits
            in_nbapp = trait_name in nbapp_traits
            if trait_name in IGNORED_TRAITS:
                # Pass trait through without any warning message.
                extapp_config_shim.update({trait_name: trait_value})
            elif all([in_extapp, in_svapp, in_nbapp]):
                warning_msg = EXTAPP_AND_NBAPP_AND_SVAPP_SHIM_MSG(
                    trait_name,
                    extapp_name
                )
                extapp_config_shim.update({trait_name: trait_value})
            elif in_extapp and in_svapp:
                warning_msg = EXTAPP_AND_SVAPP_SHIM_MSG(
                    trait_name,
                    extapp_name
                )
                extapp_config_shim.update({trait_name: trait_value})
            elif in_extapp and in_nbapp:
                warning_msg = EXTAPP_AND_NBAPP_SHIM_MSG(
                    trait_name,
                    extapp_name
                )
                extapp_config_shim.update({trait_name: trait_value})
            elif in_extapp:
                extapp_config_shim.update({trait_name: trait_value})
            elif in_svapp and in_nbapp:
                warning_msg = NOT_EXTAPP_NBAPP_AND_SVAPP_SHIM_MSG(
                    trait_name,
                    extapp_name
                )
                svapp_config_shim.update({trait_name: trait_value})
            elif in_svapp:
                warning_msg = EXTAPP_TO_SVAPP_SHIM_MSG(
                    trait_name,
                    extapp_name
                )
                svapp_config_shim.update({trait_name: trait_value})
            elif in_nbapp:
                warning_msg = EXTAPP_TO_NBAPP_SHIM_MSG(
                    trait_name,
                    extapp_name
                )
                nbapp_config_shim.update({trait_name: trait_value})
            else:
                raise TraitError("Trait, {}, not found.".format(trait_name))

            # Raise warning if one is given
            if warning_msg:
                self.log.warning(warning_msg)

        # Build config for shimmed traits.
        new_config = Config({
            'NotebookApp': nbapp_config_shim,
            'ServerApp': svapp_config_shim,
        })
        if extapp_config_shim:
            new_config.update(Config({
                self.__class__.__name__: extapp_config_shim
            }))
        # Update the full config with new values
        config_shim.update(new_config)
        return config_shim
