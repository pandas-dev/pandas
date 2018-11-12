Travis is a ci service that's well-integrated with GitHub.
The following types of breakage should be detected
by Travis builds:

1) Failing tests on any supported version of Python.
2) Pandas should install and the tests should run if no optional deps are installed.
That also means tests which rely on optional deps need to raise SkipTest()
if the dep is missing.
3) unicode related fails when running under exotic locales.

We tried running the vbench suite for a while, but with varying load
on Travis machines, that wasn't useful.

Travis currently (4/2013) has a 5-job concurrency limit. Exceeding it
basically doubles the total runtime for a commit through travis, and
since dep+pandas installation is already quite long, this should become
a hard limit on concurrent travis runs.
