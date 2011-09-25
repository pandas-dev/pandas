# reasonably effecient

def create_panels_append(cls, panels):
        """ return an appended list of panels """
        panels = [ a for a in panels if a is not None ]
        # corner cases
        if len(panels) == 0:
                return None
        elif len(panels) == 1:
                return panels[0]
        elif len(panels) == 2 and panels[0] == panels[1]:
                return panels[0]

        # add indicies that are not in the major set passed in; return a reindex version of p
        def reindex_major_axis(p, major):
                index = [ ma for ma in p.major_axis if ma not in major ]
                major.update(index)
                return p.reindex(major = index, copy = False)

        # make sure that we can actually append, e.g. that we have non-overlapping major_axis
        #   if we do, reindex so we don't
        major = set()
        panels = [ reindex_major_axis(p, major) for p in panels ]
        try:
                major =  np.concatenate([ p.major_axis for p in panels ])
        except (Exception), detail:
                raise Exception("cannot append major_axis that dont' match dimensions! -> %s" % str(detail))

        # reindex on minor axis/items
        try:
                minor, items = set(), set()
                for p in panels:
                        items.update(p.items)
                        minor.update(p.minor_axis)
                minor   = Index(sorted(list(minor)))
                items   = Index(sorted(list(items)))
                panels  = [ p.reindex(items = items, minor = minor, copy = False) for p in panels ]
        except (Exception), detail:
                raise Exception("cannot append minor/items that dont' match dimensions! -> [%s] %s" % (','.join([ "%s" % p for p in panels ]),str(detail)))

        # concatenate values
        try:
                values = np.concatenate([ p.values for p in panels ],axis=1)
        except (Exception), detail:
                raise Exception("cannot append values that dont' match dimensions! -> [%s] %s" % (','.join([ "%s" % p for p in panels ]),str(detail)))
        return Panel(values, items, major, minor )
add_class_method(Panel, create_panels_append, 'append_many')




# does the job but inefficient (better to handle like you read a table in pytables...e.g create a LongPanel then convert to Wide)

def create_panels_join(cls, panels):
        """ given an array of panels's, create a single panel """
        panels = [ a for a in panels if a is not None ]
        # corner cases
        if len(panels) == 0:
                return None
        elif len(panels) == 1:
                return panels[0]
        elif len(panels) == 2 and panels[0] == panels[1]:
                return panels[0]
        d = dict()
        minor, major, items = set(), set(), set()
        for panel in panels:
                items.update(panel.items)
                major.update(panel.major_axis)
                minor.update(panel.minor_axis)
                values = panel.values
                for item, item_index in panel.items.indexMap.items():
                        for minor_i, minor_index in panel.minor_axis.indexMap.items():
                                for major_i, major_index in panel.major_axis.indexMap.items():
                                        try:
                                                d[(minor_i,major_i,item)] = values[item_index,major_index,minor_index]
                                        except:
                                                pass
        # stack the values
        minor = sorted(list(minor))
        major = sorted(list(major))
        items = sorted(list(items))
        # create the 3d stack (items x columns x indicies)
        data = np.dstack([ np.asarray([ np.asarray([ d.get((minor_i,major_i,item),np.nan) for item in items ]) for major_i in major ]).transpose() for minor_i in minor ])
        # construct the panel
        return Panel(data, items, major, minor)
add_class_method(Panel, create_panels_join, 'join_many')

