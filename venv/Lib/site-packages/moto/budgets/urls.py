from .responses import BudgetsResponse

url_bases = [
    r"https?://budgets\.amazonaws\.com",
]


url_paths = {
    "{0}/$": BudgetsResponse.dispatch,
}
