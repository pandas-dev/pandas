# URI for the CSP Report. Included here to prevent a cyclic dependency.
# csp_report_uri is needed both by the BaseHandler (for setting the report-uri)
# and by the CSPReportHandler (which depends on the BaseHandler).
csp_report_uri = r"/api/security/csp-report"
