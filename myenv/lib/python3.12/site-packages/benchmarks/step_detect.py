try:
    from asv import step_detect
except ImportError:
    pass


class Simple:
    def setup(self):
        self.y = ([1] * 20 + [2] * 30) * 50

    if hasattr(step_detect, 'detect_steps'):
        def time_detect_regressions(self):
            steps = step_detect.detect_steps(self.y)
            step_detect.detect_regressions(steps)
    else:
        def time_detect_regressions(self):
            step_detect.detect_regressions(self.y)

    def time_solve_potts_approx(self):
        step_detect.solve_potts_approx(self.y, [1] * len(self.y), gamma=0.3)
