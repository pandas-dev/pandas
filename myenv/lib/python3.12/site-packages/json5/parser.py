# pylint: disable=line-too-long,too-many-lines,unnecessary-lambda

import unicodedata


class Parser:
    def __init__(self, msg, fname):
        self.msg = msg
        self.end = len(self.msg)
        self.fname = fname
        self.val = None
        self.pos = 0
        self.failed = False
        self.errpos = 0
        self._scopes = []
        self._cache = {}

    def parse(self):
        self._grammar_()
        if self.failed:
            return None, self._err_str(), self.errpos
        return self.val, None, self.pos

    def _err_str(self):
        lineno, colno = self._err_offsets()
        if self.errpos == len(self.msg):
            thing = 'end of input'
        else:
            thing = f'"{self.msg[self.errpos]}"'
        return f'{self.fname}:{lineno} Unexpected {thing} at column {colno}'

    def _err_offsets(self):
        lineno = 1
        colno = 1
        for i in range(self.errpos):
            if self.msg[i] == '\n':
                lineno += 1
                colno = 1
            else:
                colno += 1
        return lineno, colno

    def _succeed(self, v, newpos=None):
        self.val = v
        self.failed = False
        if newpos is not None:
            self.pos = newpos

    def _fail(self):
        self.val = None
        self.failed = True
        self.errpos = max(self.errpos, self.pos)

    def _rewind(self, newpos):
        self._succeed(None, newpos)

    def _bind(self, rule, var):
        rule()
        if not self.failed:
            self._set(var, self.val)

    def _not(self, rule):
        p = self.pos
        errpos = self.errpos
        rule()
        if self.failed:
            self._succeed(None, p)
        else:
            self._rewind(p)
            self.errpos = errpos
            self._fail()

    def _opt(self, rule):
        p = self.pos
        rule()
        if self.failed:
            self._succeed([], p)
        else:
            self._succeed([self.val])

    def _plus(self, rule):
        vs = []
        rule()
        vs.append(self.val)
        if self.failed:
            return
        self._star(rule, vs)

    def _star(self, rule, vs=None):
        vs = vs or []
        while not self.failed:
            p = self.pos
            rule()
            if self.failed:
                self._rewind(p)
                break
            vs.append(self.val)
        self._succeed(vs)

    def _seq(self, rules):
        for rule in rules:
            rule()
            if self.failed:
                return

    def _choose(self, rules):
        p = self.pos
        for rule in rules[:-1]:
            rule()
            if not self.failed:
                return
            self._rewind(p)
        rules[-1]()

    def _ch(self, ch):
        p = self.pos
        if p < self.end and self.msg[p] == ch:
            self._succeed(ch, self.pos + 1)
        else:
            self._fail()

    def _str(self, s):
        for ch in s:
            self._ch(ch)
            if self.failed:
                return
        self.val = s

    def _range(self, i, j):
        p = self.pos
        if p != self.end and ord(i) <= ord(self.msg[p]) <= ord(j):
            self._succeed(self.msg[p], self.pos + 1)
        else:
            self._fail()

    def _push(self, name):
        self._scopes.append((name, {}))

    def _pop(self, name):
        actual_name, _ = self._scopes.pop()
        assert name == actual_name

    def _get(self, var):
        return self._scopes[-1][1][var]

    def _set(self, var, val):
        self._scopes[-1][1][var] = val

    def _is_unicat(self, var, cat):
        return unicodedata.category(var) == cat

    def _join(self, s, vs):
        return s.join(vs)

    def _xtou(self, s):
        return chr(int(s, base=16))

    def _grammar_(self):
        self._push('grammar')
        self._seq(
            [
                self._sp_,
                lambda: self._bind(self._value_, 'v'),
                self._sp_,
                self._end_,
                lambda: self._succeed(self._get('v')),
            ]
        )
        self._pop('grammar')

    def _sp_(self):
        self._star(self._ws_)

    def _ws_(self):
        self._choose(
            [
                self._ws__c0_,
                self._eol_,
                self._comment_,
                self._ws__c3_,
                self._ws__c4_,
                self._ws__c5_,
                self._ws__c6_,
                self._ws__c7_,
                self._ws__c8_,
            ]
        )

    def _ws__c0_(self):
        self._ch(' ')

    def _ws__c3_(self):
        self._ch('\t')

    def _ws__c4_(self):
        self._ch('\v')

    def _ws__c5_(self):
        self._ch('\f')

    def _ws__c6_(self):
        self._ch('\xa0')

    def _ws__c7_(self):
        self._ch('\ufeff')

    def _ws__c8_(self):
        self._push('ws__c8')
        self._seq(
            [
                self._ws__c8__s0_,
                lambda: self._bind(self._anything_, 'x'),
                lambda: self._succeed(self._get('x')),
            ]
        )
        self._pop('ws__c8')

    def _ws__c8__s0_(self):
        self._not(lambda: self._not(self._ws__c8__s0_n_n_))

    def _ws__c8__s0_n_n_(self):
        self._choose([self._ws__c8__s0_n_n_g__c0_])

    def _ws__c8__s0_n_n_g__c0_(self):
        self._seq(
            [
                lambda: self._bind(self._anything_, 'x'),
                self._ws__c8__s0_n_n_g__c0__s1_,
            ]
        )

    def _ws__c8__s0_n_n_g__c0__s1_(self):
        v = self._is_unicat(self._get('x'), 'Zs')
        if v:
            self._succeed(v)
        else:
            self._fail()

    def _eol_(self):
        self._choose(
            [
                self._eol__c0_,
                self._eol__c1_,
                self._eol__c2_,
                self._eol__c3_,
                self._eol__c4_,
            ]
        )

    def _eol__c0_(self):
        self._seq([lambda: self._ch('\r'), lambda: self._ch('\n')])

    def _eol__c1_(self):
        self._ch('\r')

    def _eol__c2_(self):
        self._ch('\n')

    def _eol__c3_(self):
        self._ch('\u2028')

    def _eol__c4_(self):
        self._ch('\u2029')

    def _comment_(self):
        self._choose([self._comment__c0_, self._comment__c1_])

    def _comment__c0_(self):
        self._seq(
            [
                lambda: self._str('//'),
                lambda: self._star(self._comment__c0__s1_p_),
            ]
        )

    def _comment__c0__s1_p_(self):
        self._seq([lambda: self._not(self._eol_), self._anything_])

    def _comment__c1_(self):
        self._seq(
            [
                lambda: self._str('/*'),
                self._comment__c1__s1_,
                lambda: self._str('*/'),
            ]
        )

    def _comment__c1__s1_(self):
        self._star(
            lambda: self._seq([self._comment__c1__s1_p__s0_, self._anything_])
        )

    def _comment__c1__s1_p__s0_(self):
        self._not(lambda: self._str('*/'))

    def _value_(self):
        self._choose(
            [
                self._value__c0_,
                self._value__c1_,
                self._value__c2_,
                self._value__c3_,
                self._value__c4_,
                self._value__c5_,
                self._value__c6_,
            ]
        )

    def _value__c0_(self):
        self._seq([lambda: self._str('null'), lambda: self._succeed('None')])

    def _value__c1_(self):
        self._seq([lambda: self._str('true'), lambda: self._succeed('True')])

    def _value__c2_(self):
        self._seq([lambda: self._str('false'), lambda: self._succeed('False')])

    def _value__c3_(self):
        self._push('value__c3')
        self._seq(
            [
                lambda: self._bind(self._object_, 'v'),
                lambda: self._succeed(['object', self._get('v')]),
            ]
        )
        self._pop('value__c3')

    def _value__c4_(self):
        self._push('value__c4')
        self._seq(
            [
                lambda: self._bind(self._array_, 'v'),
                lambda: self._succeed(['array', self._get('v')]),
            ]
        )
        self._pop('value__c4')

    def _value__c5_(self):
        self._push('value__c5')
        self._seq(
            [
                lambda: self._bind(self._string_, 'v'),
                lambda: self._succeed(['string', self._get('v')]),
            ]
        )
        self._pop('value__c5')

    def _value__c6_(self):
        self._push('value__c6')
        self._seq(
            [
                lambda: self._bind(self._num_literal_, 'v'),
                lambda: self._succeed(['number', self._get('v')]),
            ]
        )
        self._pop('value__c6')

    def _object_(self):
        self._choose([self._object__c0_, self._object__c1_])

    def _object__c0_(self):
        self._push('object__c0')
        self._seq(
            [
                lambda: self._ch('{'),
                self._sp_,
                lambda: self._bind(self._member_list_, 'v'),
                self._sp_,
                lambda: self._ch('}'),
                lambda: self._succeed(self._get('v')),
            ]
        )
        self._pop('object__c0')

    def _object__c1_(self):
        self._seq(
            [
                lambda: self._ch('{'),
                self._sp_,
                lambda: self._ch('}'),
                lambda: self._succeed([]),
            ]
        )

    def _array_(self):
        self._choose([self._array__c0_, self._array__c1_])

    def _array__c0_(self):
        self._push('array__c0')
        self._seq(
            [
                lambda: self._ch('['),
                self._sp_,
                lambda: self._bind(self._element_list_, 'v'),
                self._sp_,
                lambda: self._ch(']'),
                lambda: self._succeed(self._get('v')),
            ]
        )
        self._pop('array__c0')

    def _array__c1_(self):
        self._seq(
            [
                lambda: self._ch('['),
                self._sp_,
                lambda: self._ch(']'),
                lambda: self._succeed([]),
            ]
        )

    def _string_(self):
        self._choose([self._string__c0_, self._string__c1_])

    def _string__c0_(self):
        self._push('string__c0')
        self._seq(
            [
                self._squote_,
                self._string__c0__s1_,
                self._squote_,
                lambda: self._succeed(self._join('', self._get('cs'))),
            ]
        )
        self._pop('string__c0')

    def _string__c0__s1_(self):
        self._bind(lambda: self._star(self._sqchar_), 'cs')

    def _string__c1_(self):
        self._push('string__c1')
        self._seq(
            [
                self._dquote_,
                self._string__c1__s1_,
                self._dquote_,
                lambda: self._succeed(self._join('', self._get('cs'))),
            ]
        )
        self._pop('string__c1')

    def _string__c1__s1_(self):
        self._bind(lambda: self._star(self._dqchar_), 'cs')

    def _sqchar_(self):
        self._choose([self._sqchar__c0_, self._sqchar__c1_, self._sqchar__c2_])

    def _sqchar__c0_(self):
        self._push('sqchar__c0')
        self._seq(
            [
                self._bslash_,
                lambda: self._bind(self._esc_char_, 'c'),
                lambda: self._succeed(self._get('c')),
            ]
        )
        self._pop('sqchar__c0')

    def _sqchar__c1_(self):
        self._seq([self._bslash_, self._eol_, lambda: self._succeed('')])

    def _sqchar__c2_(self):
        self._push('sqchar__c2')
        self._seq(
            [
                lambda: self._not(self._bslash_),
                lambda: self._not(self._squote_),
                lambda: self._not(self._eol_),
                lambda: self._bind(self._anything_, 'c'),
                lambda: self._succeed(self._get('c')),
            ]
        )
        self._pop('sqchar__c2')

    def _dqchar_(self):
        self._choose([self._dqchar__c0_, self._dqchar__c1_, self._dqchar__c2_])

    def _dqchar__c0_(self):
        self._push('dqchar__c0')
        self._seq(
            [
                self._bslash_,
                lambda: self._bind(self._esc_char_, 'c'),
                lambda: self._succeed(self._get('c')),
            ]
        )
        self._pop('dqchar__c0')

    def _dqchar__c1_(self):
        self._seq([self._bslash_, self._eol_, lambda: self._succeed('')])

    def _dqchar__c2_(self):
        self._push('dqchar__c2')
        self._seq(
            [
                lambda: self._not(self._bslash_),
                lambda: self._not(self._dquote_),
                lambda: self._not(self._eol_),
                lambda: self._bind(self._anything_, 'c'),
                lambda: self._succeed(self._get('c')),
            ]
        )
        self._pop('dqchar__c2')

    def _bslash_(self):
        self._ch('\\')

    def _squote_(self):
        self._ch("'")

    def _dquote_(self):
        self._ch('"')

    def _esc_char_(self):
        self._choose(
            [
                self._esc_char__c0_,
                self._esc_char__c1_,
                self._esc_char__c2_,
                self._esc_char__c3_,
                self._esc_char__c4_,
                self._esc_char__c5_,
                self._esc_char__c6_,
                self._esc_char__c7_,
                self._esc_char__c8_,
                self._esc_char__c9_,
                self._esc_char__c10_,
                self._esc_char__c11_,
                self._esc_char__c12_,
            ]
        )

    def _esc_char__c0_(self):
        self._seq([lambda: self._ch('b'), lambda: self._succeed('\b')])

    def _esc_char__c1_(self):
        self._seq([lambda: self._ch('f'), lambda: self._succeed('\f')])

    def _esc_char__c10_(self):
        self._seq(
            [
                lambda: self._ch('0'),
                lambda: self._not(self._digit_),
                lambda: self._succeed('\x00'),
            ]
        )

    def _esc_char__c11_(self):
        self._push('esc_char__c11')
        self._seq(
            [
                lambda: self._bind(self._hex_esc_, 'c'),
                lambda: self._succeed(self._get('c')),
            ]
        )
        self._pop('esc_char__c11')

    def _esc_char__c12_(self):
        self._push('esc_char__c12')
        self._seq(
            [
                lambda: self._bind(self._unicode_esc_, 'c'),
                lambda: self._succeed(self._get('c')),
            ]
        )
        self._pop('esc_char__c12')

    def _esc_char__c2_(self):
        self._seq([lambda: self._ch('n'), lambda: self._succeed('\n')])

    def _esc_char__c3_(self):
        self._seq([lambda: self._ch('r'), lambda: self._succeed('\r')])

    def _esc_char__c4_(self):
        self._seq([lambda: self._ch('t'), lambda: self._succeed('\t')])

    def _esc_char__c5_(self):
        self._seq([lambda: self._ch('v'), lambda: self._succeed('\v')])

    def _esc_char__c6_(self):
        self._seq([self._squote_, lambda: self._succeed("'")])

    def _esc_char__c7_(self):
        self._seq([self._dquote_, lambda: self._succeed('"')])

    def _esc_char__c8_(self):
        self._seq([self._bslash_, lambda: self._succeed('\\')])

    def _esc_char__c9_(self):
        self._push('esc_char__c9')
        self._seq(
            [
                self._esc_char__c9__s0_,
                lambda: self._bind(self._anything_, 'c'),
                lambda: self._succeed(self._get('c')),
            ]
        )
        self._pop('esc_char__c9')

    def _esc_char__c9__s0_(self):
        self._not(lambda: (self._esc_char__c9__s0_n_g_)())

    def _esc_char__c9__s0_n_g_(self):
        self._choose(
            [
                self._esc_char__c9__s0_n_g__c0_,
                self._esc_char__c9__s0_n_g__c1_,
                lambda: self._seq([self._digit_]),
                lambda: self._seq([self._eol_]),
            ]
        )

    def _esc_char__c9__s0_n_g__c0_(self):
        self._seq([lambda: self._ch('x')])

    def _esc_char__c9__s0_n_g__c1_(self):
        self._seq([lambda: self._ch('u')])

    def _hex_esc_(self):
        self._push('hex_esc')
        self._seq(
            [
                lambda: self._ch('x'),
                lambda: self._bind(self._hex_, 'h1'),
                lambda: self._bind(self._hex_, 'h2'),
                lambda: self._succeed(
                    self._xtou(self._get('h1') + self._get('h2'))
                ),
            ]
        )
        self._pop('hex_esc')

    def _unicode_esc_(self):
        self._push('unicode_esc')
        self._seq(
            [
                lambda: self._ch('u'),
                lambda: self._bind(self._hex_, 'a'),
                lambda: self._bind(self._hex_, 'b'),
                lambda: self._bind(self._hex_, 'c'),
                lambda: self._bind(self._hex_, 'd'),
                lambda: self._succeed(
                    self._xtou(
                        self._get('a')
                        + self._get('b')
                        + self._get('c')
                        + self._get('d')
                    )
                ),
            ]
        )
        self._pop('unicode_esc')

    def _element_list_(self):
        self._push('element_list')
        self._seq(
            [
                lambda: self._bind(self._value_, 'v'),
                self._element_list__s1_,
                self._sp_,
                self._element_list__s3_,
                lambda: self._succeed([self._get('v')] + self._get('vs')),
            ]
        )
        self._pop('element_list')

    def _element_list__s1_(self):
        self._bind(lambda: self._star(self._element_list__s1_l_p_), 'vs')

    def _element_list__s1_l_p_(self):
        self._seq([self._sp_, lambda: self._ch(','), self._sp_, self._value_])

    def _element_list__s3_(self):
        self._opt(lambda: self._ch(','))

    def _member_list_(self):
        self._push('member_list')
        self._seq(
            [
                lambda: self._bind(self._member_, 'm'),
                self._member_list__s1_,
                self._sp_,
                self._member_list__s3_,
                lambda: self._succeed([self._get('m')] + self._get('ms')),
            ]
        )
        self._pop('member_list')

    def _member_list__s1_(self):
        self._bind(lambda: self._star(self._member_list__s1_l_p_), 'ms')

    def _member_list__s1_l_p_(self):
        self._seq([self._sp_, lambda: self._ch(','), self._sp_, self._member_])

    def _member_list__s3_(self):
        self._opt(lambda: self._ch(','))

    def _member_(self):
        self._choose([self._member__c0_, self._member__c1_])

    def _member__c0_(self):
        self._push('member__c0')
        self._seq(
            [
                lambda: self._bind(self._string_, 'k'),
                self._sp_,
                lambda: self._ch(':'),
                self._sp_,
                lambda: self._bind(self._value_, 'v'),
                lambda: self._succeed([self._get('k'), self._get('v')]),
            ]
        )
        self._pop('member__c0')

    def _member__c1_(self):
        self._push('member__c1')
        self._seq(
            [
                lambda: self._bind(self._ident_, 'k'),
                self._sp_,
                lambda: self._ch(':'),
                self._sp_,
                lambda: self._bind(self._value_, 'v'),
                lambda: self._succeed([self._get('k'), self._get('v')]),
            ]
        )
        self._pop('member__c1')

    def _ident_(self):
        self._push('ident')
        self._seq(
            [
                lambda: self._bind(self._id_start_, 'hd'),
                self._ident__s1_,
                lambda: self._succeed(
                    self._join('', [self._get('hd')] + self._get('tl'))
                ),
            ]
        )
        self._pop('ident')

    def _ident__s1_(self):
        self._bind(lambda: self._star(self._id_continue_), 'tl')

    def _id_start_(self):
        self._choose(
            [self._ascii_id_start_, self._other_id_start_, self._id_start__c2_]
        )

    def _id_start__c2_(self):
        self._seq([self._bslash_, self._unicode_esc_])

    def _ascii_id_start_(self):
        self._choose(
            [
                self._ascii_id_start__c0_,
                self._ascii_id_start__c1_,
                self._ascii_id_start__c2_,
                self._ascii_id_start__c3_,
            ]
        )

    def _ascii_id_start__c0_(self):
        self._range('a', 'z')

    def _ascii_id_start__c1_(self):
        self._range('A', 'Z')

    def _ascii_id_start__c2_(self):
        self._ch('$')

    def _ascii_id_start__c3_(self):
        self._ch('_')

    def _other_id_start_(self):
        self._choose(
            [
                self._other_id_start__c0_,
                self._other_id_start__c1_,
                self._other_id_start__c2_,
                self._other_id_start__c3_,
                self._other_id_start__c4_,
                self._other_id_start__c5_,
            ]
        )

    def _other_id_start__c0_(self):
        self._push('other_id_start__c0')
        self._seq(
            [
                lambda: self._bind(self._anything_, 'x'),
                self._other_id_start__c0__s1_,
                lambda: self._succeed(self._get('x')),
            ]
        )
        self._pop('other_id_start__c0')

    def _other_id_start__c0__s1_(self):
        v = self._is_unicat(self._get('x'), 'Ll')
        if v:
            self._succeed(v)
        else:
            self._fail()

    def _other_id_start__c1_(self):
        self._push('other_id_start__c1')
        self._seq(
            [
                lambda: self._bind(self._anything_, 'x'),
                self._other_id_start__c1__s1_,
                lambda: self._succeed(self._get('x')),
            ]
        )
        self._pop('other_id_start__c1')

    def _other_id_start__c1__s1_(self):
        v = self._is_unicat(self._get('x'), 'Lm')
        if v:
            self._succeed(v)
        else:
            self._fail()

    def _other_id_start__c2_(self):
        self._push('other_id_start__c2')
        self._seq(
            [
                lambda: self._bind(self._anything_, 'x'),
                self._other_id_start__c2__s1_,
                lambda: self._succeed(self._get('x')),
            ]
        )
        self._pop('other_id_start__c2')

    def _other_id_start__c2__s1_(self):
        v = self._is_unicat(self._get('x'), 'Lo')
        if v:
            self._succeed(v)
        else:
            self._fail()

    def _other_id_start__c3_(self):
        self._push('other_id_start__c3')
        self._seq(
            [
                lambda: self._bind(self._anything_, 'x'),
                self._other_id_start__c3__s1_,
                lambda: self._succeed(self._get('x')),
            ]
        )
        self._pop('other_id_start__c3')

    def _other_id_start__c3__s1_(self):
        v = self._is_unicat(self._get('x'), 'Lt')
        if v:
            self._succeed(v)
        else:
            self._fail()

    def _other_id_start__c4_(self):
        self._push('other_id_start__c4')
        self._seq(
            [
                lambda: self._bind(self._anything_, 'x'),
                self._other_id_start__c4__s1_,
                lambda: self._succeed(self._get('x')),
            ]
        )
        self._pop('other_id_start__c4')

    def _other_id_start__c4__s1_(self):
        v = self._is_unicat(self._get('x'), 'Lu')
        if v:
            self._succeed(v)
        else:
            self._fail()

    def _other_id_start__c5_(self):
        self._push('other_id_start__c5')
        self._seq(
            [
                lambda: self._bind(self._anything_, 'x'),
                self._other_id_start__c5__s1_,
                lambda: self._succeed(self._get('x')),
            ]
        )
        self._pop('other_id_start__c5')

    def _other_id_start__c5__s1_(self):
        v = self._is_unicat(self._get('x'), 'Nl')
        if v:
            self._succeed(v)
        else:
            self._fail()

    def _id_continue_(self):
        self._choose(
            [
                self._ascii_id_start_,
                self._digit_,
                self._other_id_start_,
                self._id_continue__c3_,
                self._id_continue__c4_,
                self._id_continue__c5_,
                self._id_continue__c6_,
                self._id_continue__c7_,
                self._id_continue__c8_,
                self._id_continue__c9_,
            ]
        )

    def _id_continue__c3_(self):
        self._push('id_continue__c3')
        self._seq(
            [
                lambda: self._bind(self._anything_, 'x'),
                self._id_continue__c3__s1_,
                lambda: self._succeed(self._get('x')),
            ]
        )
        self._pop('id_continue__c3')

    def _id_continue__c3__s1_(self):
        v = self._is_unicat(self._get('x'), 'Mn')
        if v:
            self._succeed(v)
        else:
            self._fail()

    def _id_continue__c4_(self):
        self._push('id_continue__c4')
        self._seq(
            [
                lambda: self._bind(self._anything_, 'x'),
                self._id_continue__c4__s1_,
                lambda: self._succeed(self._get('x')),
            ]
        )
        self._pop('id_continue__c4')

    def _id_continue__c4__s1_(self):
        v = self._is_unicat(self._get('x'), 'Mc')
        if v:
            self._succeed(v)
        else:
            self._fail()

    def _id_continue__c5_(self):
        self._push('id_continue__c5')
        self._seq(
            [
                lambda: self._bind(self._anything_, 'x'),
                self._id_continue__c5__s1_,
                lambda: self._succeed(self._get('x')),
            ]
        )
        self._pop('id_continue__c5')

    def _id_continue__c5__s1_(self):
        v = self._is_unicat(self._get('x'), 'Nd')
        if v:
            self._succeed(v)
        else:
            self._fail()

    def _id_continue__c6_(self):
        self._push('id_continue__c6')
        self._seq(
            [
                lambda: self._bind(self._anything_, 'x'),
                self._id_continue__c6__s1_,
                lambda: self._succeed(self._get('x')),
            ]
        )
        self._pop('id_continue__c6')

    def _id_continue__c6__s1_(self):
        v = self._is_unicat(self._get('x'), 'Pc')
        if v:
            self._succeed(v)
        else:
            self._fail()

    def _id_continue__c7_(self):
        self._seq([self._bslash_, self._unicode_esc_])

    def _id_continue__c8_(self):
        self._ch('\u200c')

    def _id_continue__c9_(self):
        self._ch('\u200d')

    def _num_literal_(self):
        self._choose(
            [
                self._num_literal__c0_,
                self._num_literal__c1_,
                self._num_literal__c2_,
                self._hex_literal_,
                self._num_literal__c4_,
                self._num_literal__c5_,
            ]
        )

    def _num_literal__c0_(self):
        self._push('num_literal__c0')
        self._seq(
            [
                lambda: self._ch('-'),
                lambda: self._bind(self._num_literal_, 'n'),
                lambda: self._succeed('-' + self._get('n')),
            ]
        )
        self._pop('num_literal__c0')

    def _num_literal__c1_(self):
        self._push('num_literal__c1')
        self._seq(
            [
                lambda: self._ch('+'),
                lambda: self._bind(self._num_literal_, 'n'),
                lambda: self._succeed(self._get('n')),
            ]
        )
        self._pop('num_literal__c1')

    def _num_literal__c2_(self):
        self._push('num_literal__c2')
        self._seq(
            [
                lambda: self._bind(self._dec_literal_, 'd'),
                lambda: self._not(self._id_start_),
                lambda: self._succeed(self._get('d')),
            ]
        )
        self._pop('num_literal__c2')

    def _num_literal__c4_(self):
        self._str('Infinity')

    def _num_literal__c5_(self):
        self._str('NaN')

    def _dec_literal_(self):
        self._choose(
            [
                self._dec_literal__c0_,
                self._dec_literal__c1_,
                self._dec_literal__c2_,
                self._dec_literal__c3_,
                self._dec_literal__c4_,
                self._dec_literal__c5_,
            ]
        )

    def _dec_literal__c0_(self):
        self._push('dec_literal__c0')
        self._seq(
            [
                lambda: self._bind(self._dec_int_lit_, 'd'),
                lambda: self._bind(self._frac_, 'f'),
                lambda: self._bind(self._exp_, 'e'),
                lambda: self._succeed(
                    self._get('d') + self._get('f') + self._get('e')
                ),
            ]
        )
        self._pop('dec_literal__c0')

    def _dec_literal__c1_(self):
        self._push('dec_literal__c1')
        self._seq(
            [
                lambda: self._bind(self._dec_int_lit_, 'd'),
                lambda: self._bind(self._frac_, 'f'),
                lambda: self._succeed(self._get('d') + self._get('f')),
            ]
        )
        self._pop('dec_literal__c1')

    def _dec_literal__c2_(self):
        self._push('dec_literal__c2')
        self._seq(
            [
                lambda: self._bind(self._dec_int_lit_, 'd'),
                lambda: self._bind(self._exp_, 'e'),
                lambda: self._succeed(self._get('d') + self._get('e')),
            ]
        )
        self._pop('dec_literal__c2')

    def _dec_literal__c3_(self):
        self._push('dec_literal__c3')
        self._seq(
            [
                lambda: self._bind(self._dec_int_lit_, 'd'),
                lambda: self._succeed(self._get('d')),
            ]
        )
        self._pop('dec_literal__c3')

    def _dec_literal__c4_(self):
        self._push('dec_literal__c4')
        self._seq(
            [
                lambda: self._bind(self._frac_, 'f'),
                lambda: self._bind(self._exp_, 'e'),
                lambda: self._succeed(self._get('f') + self._get('e')),
            ]
        )
        self._pop('dec_literal__c4')

    def _dec_literal__c5_(self):
        self._push('dec_literal__c5')
        self._seq(
            [
                lambda: self._bind(self._frac_, 'f'),
                lambda: self._succeed(self._get('f')),
            ]
        )
        self._pop('dec_literal__c5')

    def _dec_int_lit_(self):
        self._choose([self._dec_int_lit__c0_, self._dec_int_lit__c1_])

    def _dec_int_lit__c0_(self):
        self._seq(
            [
                lambda: self._ch('0'),
                lambda: self._not(self._digit_),
                lambda: self._succeed('0'),
            ]
        )

    def _dec_int_lit__c1_(self):
        self._push('dec_int_lit__c1')
        self._seq(
            [
                lambda: self._bind(self._nonzerodigit_, 'd'),
                self._dec_int_lit__c1__s1_,
                lambda: self._succeed(
                    self._get('d') + self._join('', self._get('ds'))
                ),
            ]
        )
        self._pop('dec_int_lit__c1')

    def _dec_int_lit__c1__s1_(self):
        self._bind(lambda: self._star(self._digit_), 'ds')

    def _digit_(self):
        self._range('0', '9')

    def _nonzerodigit_(self):
        self._range('1', '9')

    def _hex_literal_(self):
        self._push('hex_literal')
        self._seq(
            [
                self._hex_literal__s0_,
                self._hex_literal__s1_,
                lambda: self._succeed('0x' + self._join('', self._get('hs'))),
            ]
        )
        self._pop('hex_literal')

    def _hex_literal__s0_(self):
        self._choose([lambda: self._str('0x'), lambda: self._str('0X')])

    def _hex_literal__s1_(self):
        self._bind(lambda: self._plus(self._hex_), 'hs')

    def _hex_(self):
        self._choose([self._hex__c0_, self._hex__c1_, self._digit_])

    def _hex__c0_(self):
        self._range('a', 'f')

    def _hex__c1_(self):
        self._range('A', 'F')

    def _frac_(self):
        self._push('frac')
        self._seq(
            [
                lambda: self._ch('.'),
                self._frac__s1_,
                lambda: self._succeed('.' + self._join('', self._get('ds'))),
            ]
        )
        self._pop('frac')

    def _frac__s1_(self):
        self._bind(lambda: self._star(self._digit_), 'ds')

    def _exp_(self):
        self._choose([self._exp__c0_, self._exp__c1_])

    def _exp__c0_(self):
        self._push('exp__c0')
        self._seq(
            [
                self._exp__c0__s0_,
                lambda: self._bind(self._exp__c0__s1_l_, 's'),
                self._exp__c0__s2_,
                lambda: self._succeed(
                    'e' + self._get('s') + self._join('', self._get('ds'))
                ),
            ]
        )
        self._pop('exp__c0')

    def _exp__c0__s0_(self):
        self._choose([lambda: self._ch('e'), lambda: self._ch('E')])

    def _exp__c0__s1_l_(self):
        self._choose([lambda: self._ch('+'), lambda: self._ch('-')])

    def _exp__c0__s2_(self):
        self._bind(lambda: self._star(self._digit_), 'ds')

    def _exp__c1_(self):
        self._push('exp__c1')
        self._seq(
            [
                self._exp__c1__s0_,
                self._exp__c1__s1_,
                lambda: self._succeed('e' + self._join('', self._get('ds'))),
            ]
        )
        self._pop('exp__c1')

    def _exp__c1__s0_(self):
        self._choose([lambda: self._ch('e'), lambda: self._ch('E')])

    def _exp__c1__s1_(self):
        self._bind(lambda: self._star(self._digit_), 'ds')

    def _anything_(self):
        if self.pos < self.end:
            self._succeed(self.msg[self.pos], self.pos + 1)
        else:
            self._fail()

    def _end_(self):
        if self.pos == self.end:
            self._succeed(None)
        else:
            self._fail()
