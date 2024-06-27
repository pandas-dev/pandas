import json
import re
from typing import Any, Dict, Tuple, Union
from urllib.parse import urlsplit

from moto.core.responses import BaseResponse

from .models import PollyBackend, polly_backends
from .resources import LANGUAGE_CODES, VOICE_IDS

LEXICON_NAME_REGEX = re.compile(r"^[0-9A-Za-z]{1,20}$")


class PollyResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="polly")

    @property
    def polly_backend(self) -> PollyBackend:
        return polly_backends[self.current_account][self.region]

    @property
    def json(self) -> Dict[str, Any]:  # type: ignore[misc]
        if not hasattr(self, "_json"):
            self._json = json.loads(self.body)
        return self._json

    def _error(self, code: str, message: str) -> Tuple[str, Dict[str, int]]:
        return json.dumps({"__type": code, "message": message}), dict(status=400)

    def _get_action(self) -> str:
        # Amazon is now naming things /v1/api_name
        url_parts = urlsplit(self.uri).path.lstrip("/").split("/")
        # [0] = 'v1'

        return url_parts[1]

    # DescribeVoices
    def voices(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        language_code = self._get_param("LanguageCode")

        if language_code is not None and language_code not in LANGUAGE_CODES:
            all_codes = ", ".join(LANGUAGE_CODES)  # type: ignore
            msg = (
                f"1 validation error detected: Value '{language_code}' at 'languageCode' failed to satisfy constraint: "
                f"Member must satisfy enum value set: [{all_codes}]"
            )
            return msg, dict(status=400)

        voices = self.polly_backend.describe_voices(language_code)

        return json.dumps({"Voices": voices})

    def lexicons(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        # Dish out requests based on methods

        # anything after the /v1/lexicons/
        args = urlsplit(self.uri).path.lstrip("/").split("/")[2:]

        if self.method == "GET":
            if len(args) == 0:
                return self._get_lexicons_list()
            else:
                return self._get_lexicon(*args)
        elif self.method == "PUT":
            return self._put_lexicons(*args)
        elif self.method == "DELETE":
            return self._delete_lexicon(*args)

        return self._error("InvalidAction", "Bad route")

    # PutLexicon
    def _put_lexicons(
        self, lexicon_name: str
    ) -> Union[str, Tuple[str, Dict[str, int]]]:
        if LEXICON_NAME_REGEX.match(lexicon_name) is None:
            return self._error(
                "InvalidParameterValue", "Lexicon name must match [0-9A-Za-z]{1,20}"
            )

        if "Content" not in self.json:
            return self._error("MissingParameter", "Content is missing from the body")

        self.polly_backend.put_lexicon(lexicon_name, self.json["Content"])

        return ""

    # ListLexicons
    def _get_lexicons_list(self) -> str:
        result = {"Lexicons": self.polly_backend.list_lexicons()}

        return json.dumps(result)

    # GetLexicon
    def _get_lexicon(self, lexicon_name: str) -> Union[str, Tuple[str, Dict[str, int]]]:
        try:
            lexicon = self.polly_backend.get_lexicon(lexicon_name)
        except KeyError:
            return self._error("LexiconNotFoundException", "Lexicon not found")

        result = {
            "Lexicon": {"Name": lexicon_name, "Content": lexicon.content},
            "LexiconAttributes": lexicon.to_dict()["Attributes"],
        }

        return json.dumps(result)

    # DeleteLexicon
    def _delete_lexicon(
        self, lexicon_name: str
    ) -> Union[str, Tuple[str, Dict[str, int]]]:
        try:
            self.polly_backend.delete_lexicon(lexicon_name)
        except KeyError:
            return self._error("LexiconNotFoundException", "Lexicon not found")

        return ""

    # SynthesizeSpeech
    def speech(self) -> Tuple[str, Dict[str, Any]]:
        # Sanity check params
        args = {
            "lexicon_names": None,
            "sample_rate": 22050,
            "speech_marks": None,
            "text": None,
            "text_type": "text",
        }

        if "LexiconNames" in self.json:
            for lex in self.json["LexiconNames"]:
                try:
                    self.polly_backend.get_lexicon(lex)
                except KeyError:
                    return self._error("LexiconNotFoundException", "Lexicon not found")

            args["lexicon_names"] = self.json["LexiconNames"]

        if "OutputFormat" not in self.json:
            return self._error("MissingParameter", "Missing parameter OutputFormat")
        if self.json["OutputFormat"] not in ("json", "mp3", "ogg_vorbis", "pcm"):
            return self._error(
                "InvalidParameterValue", "Not one of json, mp3, ogg_vorbis, pcm"
            )
        args["output_format"] = self.json["OutputFormat"]

        if "SampleRate" in self.json:
            sample_rate = int(self.json["SampleRate"])
            if sample_rate not in (8000, 16000, 22050):
                return self._error(
                    "InvalidSampleRateException",
                    "The specified sample rate is not valid.",
                )
            args["sample_rate"] = sample_rate

        if "SpeechMarkTypes" in self.json:
            for value in self.json["SpeechMarkTypes"]:
                if value not in ("sentance", "ssml", "viseme", "word"):
                    return self._error(
                        "InvalidParameterValue",
                        "Not one of sentance, ssml, viseme, word",
                    )
            args["speech_marks"] = self.json["SpeechMarkTypes"]

        if "Text" not in self.json:
            return self._error("MissingParameter", "Missing parameter Text")
        args["text"] = self.json["Text"]

        if "TextType" in self.json:
            if self.json["TextType"] not in ("ssml", "text"):
                return self._error("InvalidParameterValue", "Not one of ssml, text")
            args["text_type"] = self.json["TextType"]

        if "VoiceId" not in self.json:
            return self._error("MissingParameter", "Missing parameter VoiceId")
        if self.json["VoiceId"] not in VOICE_IDS:
            all_voices = ", ".join(VOICE_IDS)  # type: ignore
            return self._error("InvalidParameterValue", f"Not one of {all_voices}")
        args["voice_id"] = self.json["VoiceId"]

        # More validation
        if len(args["text"]) > 3000:  # type: ignore
            return self._error("TextLengthExceededException", "Text too long")

        if args["speech_marks"] is not None and args["output_format"] != "json":
            return self._error(
                "MarksNotSupportedForFormatException", "OutputFormat must be json"
            )
        if args["speech_marks"] is not None and args["text_type"] == "text":
            return self._error(
                "SsmlMarksNotSupportedForTextTypeException", "TextType must be ssml"
            )

        content_type = "audio/json"
        if args["output_format"] == "mp3":
            content_type = "audio/mpeg"
        elif args["output_format"] == "ogg_vorbis":
            content_type = "audio/ogg"
        elif args["output_format"] == "pcm":
            content_type = "audio/pcm"

        headers = {"Content-Type": content_type}

        return "\x00\x00\x00\x00\x00\x00\x00\x00", headers
