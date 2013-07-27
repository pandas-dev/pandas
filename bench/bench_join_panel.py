# reasonably efficient


def create_panels_append(cls, panels):
        """ return an append list of panels """
        panels = [a for a in panels if a is not None]
        # corner cases
        if len(panels) == 0:
                return None
        elif len(panels) == 1:
                return panels[0]
        elif len(panels) == 2 and panels[0] == panels[1]:
                return panels[0]
        # import pdb; pdb.set_trace()
        # create a joint index for the axis

        def joint_index_for_axis(panels, axis):
                s = set()
                for p in panels:
                        s.update(list(getattr(p, axis)))
                return sorted(list(s))

        def reindex_on_axis(panels, axis, axis_reindex):
                new_axis = joint_index_for_axis(panels, axis)
                new_panels = [p.reindex(**{axis_reindex: new_axis,
                                        'copy': False}) for p in panels]
                return new_panels, new_axis
        # create the joint major index, dont' reindex the sub-panels - we are
        # appending
        major = joint_index_for_axis(panels, 'major_axis')
        # reindex on minor axis
        panels, minor = reindex_on_axis(panels, 'minor_axis', 'minor')
        # reindex on items
        panels, items = reindex_on_axis(panels, 'items', 'items')
        # concatenate values
        try:
                values = np.concatenate([p.values for p in panels], axis=1)
        except Exception as detail:
                raise Exception("cannot append values that dont' match dimensions! -> [%s] %s"
                                % (','.join(["%s" % p for p in panels]), str(detail)))
        # pm('append - create_panel')
        p = Panel(values, items=items, major_axis=major,
                  minor_axis=minor)
        # pm('append - done')
        return p


# does the job but inefficient (better to handle like you read a table in
# pytables...e.g create a LongPanel then convert to Wide)
def create_panels_join(cls, panels):
        """ given an array of panels's, create a single panel """
        panels = [a for a in panels if a is not None]
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
                                                d[(minor_i, major_i, item)] = values[item_index, major_index, minor_index]
                                        except:
                                                pass
        # stack the values
        minor = sorted(list(minor))
        major = sorted(list(major))
        items = sorted(list(items))
        # create the 3d stack (items x columns x indicies)
        data = np.dstack([np.asarray([np.asarray([d.get((minor_i, major_i, item), np.nan)
                                                  for item in items])
                                      for major_i in major]).transpose()
                          for minor_i in minor])
        # construct the panel
        return Panel(data, items, major, minor)
add_class_method(Panel, create_panels_join, 'join_many')
