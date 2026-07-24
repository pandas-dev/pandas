import tensorflow as tf
from tensorflow._aliases import Integer, StringTensorCompatible

def decode_wav(
    contents: StringTensorCompatible, desired_channels: int = -1, desired_samples: int = -1, name: str | None = None
) -> tuple[tf.Tensor, tf.Tensor]: ...
def encode_wav(audio: tf.Tensor, sample_rate: Integer, name: str | None = None) -> tf.Tensor: ...
