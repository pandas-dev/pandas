from  abc import abstractmethod
import os
import uuid


ROW_HEADING_CLASS = "pandas_row_heading"
COL_HEADING_CLASS = "pandas_col_heading"
DATA_CLASS = "pandas_data"
DATA_CELL_TYPE = "data"
ROW_CLASS = "pandas_row"
COLUMN_CLASS = "pandas_col"
HEADING_CELL_TYPE = "heading"
BLANK_CLASS = "pandas_blank"
BLANK_VALUE = ""
LEVEL_CLASS = "pandas_level"

TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "templates")

# TODO, cleaner handling of cell_context (extra_contaxt
class Styler(object):

    @abstractmethod
    def render(self, f):
        pass

    def __init__(self, df, template, engine_instance=None,*args, **kwds):
        super(Styler, self).__init__(*args, **kwds)
        self.df = df
        self.template = template
        self.style = []
        self.cell_context = {}

        from pandas.io.templating.engines import Jinja2Engine
        self.engine_instance = engine_instance or Jinja2Engine()
        self.engine_instance.load(self.template)

    def std_context(self, *args, **kwds):
        df = self.df

        n_rlvls = df.index.nlevels
        n_clvls = df.columns.nlevels
        rlabels = df.index.tolist()
        clabels = df.columns.tolist()
        if n_rlvls == 1:
            rlabels = [[x] for x in rlabels]
        if n_clvls == 1:
            clabels = [[x] for x in clabels]
        clabels = zip(*clabels)
        head = []
        for r in range(n_clvls):
            row_es = [{"type": HEADING_CELL_TYPE, "value": BLANK_VALUE,
                       "class": " ".join([BLANK_CLASS])}] * n_rlvls
            for c in range(len(clabels[0])):
                cs = [COL_HEADING_CLASS, "%s%s" % (LEVEL_CLASS, r),
                      "%s%s" % (COLUMN_CLASS, c)]
                cs.extend(
                    self.cell_context.get("col_headings", {}).get(r, {}).get(c,
                        []))
                row_es.append(
                    {"type": HEADING_CELL_TYPE, "value": clabels[r][c],
                     "class": " ".join(cs)})
            head.append(row_es)
        body = []
        for r in range(len(df)):
            cs = [ROW_HEADING_CLASS, "%s%s" % (LEVEL_CLASS, c),
                  "%s%s" % (ROW_CLASS, r)]
            cs.extend(
                self.cell_context.get("row_headings", {}).get(r, {}).get(c, []))
            row_es = [
                {"type": HEADING_CELL_TYPE, "value": rlabels[r][c],
                 "class": " ".join(cs)}
                for c in range(len(rlabels[r]))]
            for c in range(len(df.columns)):
                cs = [DATA_CLASS, "%s%s" % (ROW_CLASS, r),
                      "%s%s" % (COLUMN_CLASS, c)]
                cs.extend(
                    self.cell_context.get("data", {}).get(r, {}).get(c, []))
                row_es.append(
                    {"type": DATA_CELL_TYPE, "value": df.iloc[r][c],
                     "class": " ".join(cs)})
            body.append(row_es)

        # uuid required to isolate table styling from other tables
        # on the page in ipnb
        u = str(uuid.uuid1()).replace("-", "_")
        return dict(head=head, body=body, uuid=u, style=self.style)

    def render(self, f=None, **kwds):
        encoding = kwds.pop('encoding', "utf8")
        s = self.engine_instance.render(self.std_context())
        if f:
            with codecs.open(f, "wb", encoding) as f:
                f.write(s)
        else:
            return s

class ITemplateEngine(object):
    """Interface for supporting multiple template engines

    We'll support only a single engine, but help users help themselves.
    """


    @abstractmethod
    def render(self, f=None,**kwds):
        pass

    @abstractmethod
    def load(self, *args, **kwds):
        pass


from html import HTMLStyler