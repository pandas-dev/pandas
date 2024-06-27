from moto.moto_api import _internal

"""
Global StateManager that everyone uses
Use this manager to configure how AWS models transition between states. (initializing -> starting, starting -> ready, etc.)
"""
state_manager = _internal.state_manager.StateManager()

"""
Default transitions across Moto
"""
state_manager.register_default_transition(
    "batch::job", transition={"progression": "manual", "times": 1}
)
state_manager.register_default_transition(
    "cloudfront::distribution", transition={"progression": "manual", "times": 1}
)
state_manager.register_default_transition(
    model_name="dax::cluster", transition={"progression": "manual", "times": 4}
)
state_manager.register_default_transition(
    model_name="ecs::task", transition={"progression": "manual", "times": 1}
)
state_manager.register_default_transition(
    model_name="glue::job_run", transition={"progression": "immediate"}
)
state_manager.register_default_transition(
    "s3::keyrestore", transition={"progression": "immediate"}
)
state_manager.register_default_transition(
    model_name="support::case", transition={"progression": "manual", "times": 1}
)
state_manager.register_default_transition(
    "transcribe::vocabulary", transition={"progression": "manual", "times": 1}
)
state_manager.register_default_transition(
    "transcribe::medicalvocabulary", transition={"progression": "manual", "times": 1}
)
state_manager.register_default_transition(
    "transcribe::transcriptionjob", transition={"progression": "manual", "times": 1}
)
state_manager.register_default_transition(
    "transcribe::medicaltranscriptionjob",
    transition={"progression": "manual", "times": 1},
)

""""
Recorder, used to record calls to Moto and replay them later
"""
recorder = _internal.Recorder()
