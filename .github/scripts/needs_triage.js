// Apply `Needs Triage` to an issue opened by someone other than an owner,
// member, collaborator, or bot, unless the label has been applied or removed
// before. Returns whether the issue carries the label afterwards.
const LABEL = 'Needs Triage';
const EXEMPT_ASSOCIATIONS = ['OWNER', 'MEMBER', 'COLLABORATOR'];
// Only apply going forward; not retroactively.
const ROLLOUT_CUTOFF = new Date('????-??-??T00:00:00Z');

module.exports = async ({ github, context, issue }) => {
  const labels = issue.labels.map(label => typeof label === 'string' ? label : label.name);
  if (labels.includes(LABEL)) return true;
  if (new Date(issue.created_at) < ROLLOUT_CUTOFF) return false;
  if (issue.user.type === 'Bot' || EXEMPT_ASSOCIATIONS.includes(issue.author_association)) {
    return false;
  }
  const repo = { owner: context.repo.owner, repo: context.repo.repo };
  const events = await github.paginate(github.rest.issues.listEvents, {
    ...repo, issue_number: issue.number, per_page: 100,
  });
  const seen = events.some(event =>
    ['labeled', 'unlabeled'].includes(event.event) && event.label?.name === LABEL);
  if (seen) return false;
  await github.rest.issues.addLabels({ ...repo, issue_number: issue.number, labels: [LABEL] });
  return true;
};
