class FailBenchmarks:
    def time_cause_asv_to_exit_with_non_zero_code(self):
        raise ValueError('If the CI build is red, things work as expected')
