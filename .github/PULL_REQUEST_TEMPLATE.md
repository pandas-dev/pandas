Checklist for the pandas documentation sprint (ignore this if you are doing
an unrelated PR):

- [ ] PR title is "DOC: update the <your-function-or-method> docstring"
- [ ] The validation script passes: `scripts/validate_docstrings.py <your-function-or-method>`
- [ ] The PEP8 style check passes: `git diff upstream/master -u -- "*.py" | flake8 --diff`
- [ ] The html version looks good: `python doc/make.py --single <your-function-or-method>`
- [ ] It has been proofread on language by another sprint participant

Please include the output of the validation script below between the "```" ticks:

```
# paste output of "scripts/validate_docstrings.py <your-function-or-method>" here
# between the "```" (remove this comment, but keep the "```")

```

If the validation script still gives errors, but you think there is a good reason
to deviate in this case (and there are certainly such cases), please state this
explicitly.


Checklist for other PRs (remove this part if you are doing a PR for the pandas documentation sprint):

- [ ] closes #xxxx
- [ ] tests added / passed
- [ ] passes `git diff upstream/master -u -- "*.py" | flake8 --diff`
- [ ] whatsnew entry
