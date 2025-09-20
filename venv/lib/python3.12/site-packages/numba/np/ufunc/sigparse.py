import tokenize
import string


def parse_signature(sig):
    '''Parse generalized ufunc signature.

    NOTE: ',' (COMMA) is a delimiter; not separator.
          This means trailing comma is legal.
    '''
    def stripws(s):
        return ''.join(c for c in s if c not in string.whitespace)

    def tokenizer(src):
        def readline():
            yield src
        gen = readline()
        return tokenize.generate_tokens(lambda: next(gen))

    def parse(src):
        tokgen = tokenizer(src)
        while True:
            tok = next(tokgen)
            if tok[1] == '(':
                symbols = []
                while True:
                    tok = next(tokgen)
                    if tok[1] == ')':
                        break
                    elif tok[0] == tokenize.NAME:
                        symbols.append(tok[1])
                    elif tok[1] == ',':
                        continue
                    else:
                        raise ValueError('bad token in signature "%s"' % tok[1])
                yield tuple(symbols)
                tok = next(tokgen)
                if tok[1] == ',':
                    continue
                elif tokenize.ISEOF(tok[0]):
                    break
            elif tokenize.ISEOF(tok[0]):
                break
            else:
                raise ValueError('bad token in signature "%s"' % tok[1])

    ins, _, outs = stripws(sig).partition('->')
    inputs = list(parse(ins))
    outputs = list(parse(outs))

    # check that all output symbols are defined in the inputs
    isym = set()
    osym = set()
    for grp in inputs:
        isym |= set(grp)
    for grp in outputs:
        osym |= set(grp)

    diff = osym.difference(isym)
    if diff:
        raise NameError('undefined output symbols: %s' % ','.join(sorted(diff)))

    return inputs, outputs
