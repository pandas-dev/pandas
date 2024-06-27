#############################################################################
##
## Copyright (C) 2014 Riverbank Computing Limited.
## Copyright (C) 2006 Thorsten Marek.
## All right reserved.
##
## This file is part of PyQt.
##
## You may use this file under the terms of the GPL v2 or the revised BSD
## license as follows:
##
## "Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##   * Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##   * Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##   * Neither the name of the Riverbank Computing Limited nor the names
##     of its contributors may be used to endorse or promote products
##     derived from this software without specific prior written
##     permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
##
#############################################################################


from .misc import Literal, moduleMember


class ProxyMetaclass(type):
    """ ProxyMetaclass is the meta-class for proxies. """

    def __init__(*args):
        """ Initialise the meta-class. """

        # Initialise as normal.
        type.__init__(*args)

        # The proxy type object we have created.
        proxy = args[0]

        # Go through the proxy's attributes looking for other proxies.
        for sub_proxy in proxy.__dict__.values():
            if type(sub_proxy) is ProxyMetaclass:
                # Set the module name of the contained proxy to the name of the
                # container proxy.
                sub_proxy.module = proxy.__name__

                # Attribute hierachies are created depth first so any proxies
                # contained in the sub-proxy whose module we have just set will
                # already exist and have an incomplete module name.  We need to
                # revisit them and prepend the new name to their module names.
                # Note that this should be recursive but with current usage we
                # know there will be only one level to revisit.
                for sub_sub_proxy in sub_proxy.__dict__.values():
                    if type(sub_sub_proxy) is ProxyMetaclass:
                        sub_sub_proxy.module = '%s.%s' % (proxy.__name__, sub_sub_proxy.module)

        # Makes sure there is a 'module' attribute.
        if not hasattr(proxy, 'module'):
            proxy.module = ''
    
    def __getattribute__(cls, name):
        try:
            return type.__getattribute__(cls, name)
        except AttributeError:
            # Make sure __init__()'s use of hasattr() works.
            if name == 'module':
                raise

            # Avoid a circular import.
            from .qtproxies import LiteralProxyClass

            return type(name, (LiteralProxyClass, ),
                        {"module": moduleMember(type.__getattribute__(cls, "module"),
                                                type.__getattribute__(cls, "__name__"))})            

    def __str__(cls):
        return moduleMember(type.__getattribute__(cls, "module"),
                            type.__getattribute__(cls, "__name__"))

    def __or__(self, r_op):
        return Literal("%s|%s" % (self, r_op))

    def __eq__(self, other):
        return str(self) == str(other)
