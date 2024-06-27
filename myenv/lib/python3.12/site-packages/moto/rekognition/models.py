"""RekognitionBackend class with methods for supported APIs."""

import string
from typing import Any, Dict, List, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.moto_api._internal import mock_random as random


class RekognitionBackend(BaseBackend):
    """Implementation of Rekognition APIs."""

    def start_face_search(self) -> str:
        return self._job_id()

    def start_text_detection(self) -> str:
        return self._job_id()

    def get_face_search(
        self,
    ) -> Tuple[str, str, Dict[str, Any], List[Dict[str, Any]], str, str]:
        """
        This returns hardcoded values and none of the parameters are taken into account.
        """
        return (
            self._job_status(),
            self._status_message(),
            self._video_metadata(),
            self._persons(),
            self._next_token(),
            self._text_model_version(),
        )

    def get_text_detection(
        self,
    ) -> Tuple[str, str, Dict[str, Any], List[Dict[str, Any]], str, str]:
        """
        This returns hardcoded values and none of the parameters are taken into account.
        """
        return (
            self._job_status(),
            self._status_message(),
            self._video_metadata(),
            self._text_detections(),
            self._next_token(),
            self._text_model_version(),
        )

    def compare_faces(
        self,
    ) -> Tuple[
        List[Dict[str, Any]],
        str,
        str,
        List[Dict[str, Any]],
        Dict[str, Any],
    ]:
        return (
            self._face_matches(),
            "ROTATE_90",
            "ROTATE_90",
            self._unmatched_faces(),
            self.source_image_face(),
        )

    def detect_labels(self) -> Tuple[List[Dict[str, Any]], Dict[str, Any], str]:
        return (
            self._mobile_phone_label(),
            self._image_properties(),
            "3.0",
        )

    def detect_text(self) -> Tuple[List[Dict[str, Any]], str]:
        return (
            self._detect_text_text_detections(),
            "3.0",
        )

    def detect_custom_labels(self) -> Tuple[List[Dict[str, Any]]]:
        return (self._detect_custom_labels_detections(),)

    # private

    def _job_id(self) -> str:
        return "".join(
            random.choice(string.ascii_uppercase + string.digits) for _ in range(64)
        )

    def _job_status(self) -> str:
        return "SUCCEEDED"

    def _next_token(self) -> str:
        return ""

    def _status_message(self) -> str:
        return ""

    def _text_model_version(self) -> str:
        return "3.1"

    def _video_metadata(self) -> Dict[str, Any]:
        return {
            "Codec": "h264",
            "DurationMillis": 15020,
            "Format": "QuickTime / MOV",
            "FrameRate": 24.0,
            "FrameHeight": 720,
            "FrameWidth": 1280,
            "ColorRange": "LIMITED",
        }

    def _persons(self) -> List[Dict[str, Any]]:
        return [
            {
                "Timestamp": 0,
                "Person": {
                    "Index": 0,
                    "Face": {
                        "BoundingBox": {
                            "Width": 0.42217350006103516,
                            "Height": 0.9352386593818665,
                            "Left": 0.31870967149734497,
                            "Top": -0.0049947104416787624,
                        },
                        "Landmarks": [
                            {
                                "Type": "eyeLeft",
                                "X": 0.4800040125846863,
                                "Y": 0.23425640165805817,
                            },
                            {
                                "Type": "eyeRight",
                                "X": 0.63795405626297,
                                "Y": 0.19219470024108887,
                            },
                            {
                                "Type": "mouthLeft",
                                "X": 0.5283276438713074,
                                "Y": 0.6190487146377563,
                            },
                            {
                                "Type": "mouthRight",
                                "X": 0.660395085811615,
                                "Y": 0.5830448269844055,
                            },
                            {
                                "Type": "nose",
                                "X": 0.619724690914154,
                                "Y": 0.3800361752510071,
                            },
                        ],
                        "Pose": {
                            "Roll": -5.063229084014893,
                            "Yaw": 18.038856506347656,
                            "Pitch": 12.567241668701172,
                        },
                        "Quality": {
                            "Brightness": 83.42264556884766,
                            "Sharpness": 67.22731018066406,
                        },
                        "Confidence": 99.99860382080078,
                    },
                },
                "FaceMatches": [
                    {
                        "Similarity": 99.99994659423828,
                        "Face": {
                            "FaceId": "f2489050-020e-4c14-8693-63339847a59d",
                            "BoundingBox": {
                                "Width": 0.7136539816856384,
                                "Height": 0.9471719861030579,
                                "Left": 0.19036999344825745,
                                "Top": -0.012074699625372887,
                            },
                            "ImageId": "f3b180d3-f5ad-39c1-b825-ba30b170a90d",
                            "ExternalImageId": "Dave_Bloggs",
                            "Confidence": 99.99970245361328,
                        },
                    },
                    {
                        "Similarity": 99.9986572265625,
                        "Face": {
                            "FaceId": "f0d22a6a-3436-4d23-ae5b-c5cb2e795581",
                            "BoundingBox": {
                                "Width": 0.7198730111122131,
                                "Height": 1.003640055656433,
                                "Left": 0.1844159960746765,
                                "Top": -0.00142729002982378,
                            },
                            "ImageId": "738d14f3-26be-3066-b1a9-7f4f6bb3ffc6",
                            "ExternalImageId": "Dave_Bloggs",
                            "Confidence": 99.99939727783203,
                        },
                    },
                    {
                        "Similarity": 99.99791717529297,
                        "Face": {
                            "FaceId": "c48162bd-a16a-4e04-ad3c-967761895295",
                            "BoundingBox": {
                                "Width": 0.7364680171012878,
                                "Height": 1.0104399919509888,
                                "Left": 0.1361449956893921,
                                "Top": -0.009593159891664982,
                            },
                            "ImageId": "eae3565c-741b-342c-8e73-379a09ae5346",
                            "ExternalImageId": "Dave_Bloggs",
                            "Confidence": 99.99949645996094,
                        },
                    },
                    {
                        "Similarity": 99.37212371826172,
                        "Face": {
                            "FaceId": "651314bb-28d4-405d-9b13-c32e9ff28299",
                            "BoundingBox": {
                                "Width": 0.3711090087890625,
                                "Height": 0.3609749972820282,
                                "Left": 0.2571589946746826,
                                "Top": 0.21493400633335114,
                            },
                            "ImageId": "068700f5-0b2e-39c0-874b-2c58fa10d833",
                            "ExternalImageId": "Dave_Bloggs",
                            "Confidence": 99.99300384521484,
                        },
                    },
                ],
            }
        ]

    def _text_detections(self) -> List[Dict[str, Any]]:
        return [
            {
                "Timestamp": 0,
                "TextDetection": {
                    "DetectedText": "Hello world",
                    "Type": "LINE",
                    "Id": 0,
                    "Confidence": 97.89398956298828,
                    "Geometry": {
                        "BoundingBox": {
                            "Width": 0.1364741027355194,
                            "Height": 0.0318513885140419,
                            "Left": 0.4310702085494995,
                            "Top": 0.876121461391449,
                        },
                        "Polygon": [
                            {"X": 0.4310702085494995, "Y": 0.8769540190696716},
                            {"X": 0.5673548579216003, "Y": 0.876121461391449},
                            {"X": 0.5675443410873413, "Y": 0.90714031457901},
                            {"X": 0.4312596917152405, "Y": 0.9079728722572327},
                        ],
                    },
                },
            },
            {
                "Timestamp": 0,
                "TextDetection": {
                    "DetectedText": "Hello",
                    "Type": "WORD",
                    "Id": 1,
                    "ParentId": 0,
                    "Confidence": 99.1568832397461,
                    "Geometry": {
                        "BoundingBox": {
                            "Width": 0.0648193359375,
                            "Height": 0.0234375,
                            "Left": 0.43121337890625,
                            "Top": 0.876953125,
                        },
                        "Polygon": [
                            {"X": 0.43121337890625, "Y": 0.876953125},
                            {"X": 0.49603271484375, "Y": 0.876953125},
                            {"X": 0.49603271484375, "Y": 0.900390625},
                            {"X": 0.43121337890625, "Y": 0.900390625},
                        ],
                    },
                },
            },
            {
                "Timestamp": 0,
                "TextDetection": {
                    "DetectedText": "world",
                    "Type": "WORD",
                    "Id": 2,
                    "ParentId": 0,
                    "Confidence": 96.63108825683594,
                    "Geometry": {
                        "BoundingBox": {
                            "Width": 0.07103776931762695,
                            "Height": 0.02804870530962944,
                            "Left": 0.4965003430843353,
                            "Top": 0.8795245885848999,
                        },
                        "Polygon": [
                            {"X": 0.4965003430843353, "Y": 0.8809727430343628},
                            {"X": 0.5673661231994629, "Y": 0.8795245885848999},
                            {"X": 0.5675381422042847, "Y": 0.9061251282691956},
                            {"X": 0.4966723322868347, "Y": 0.9075732827186584},
                        ],
                    },
                },
            },
            {
                "Timestamp": 1000,
                "TextDetection": {
                    "DetectedText": "Goodbye world",
                    "Type": "LINE",
                    "Id": 0,
                    "Confidence": 98.9729995727539,
                    "Geometry": {
                        "BoundingBox": {
                            "Width": 0.13677978515625,
                            "Height": 0.0302734375,
                            "Left": 0.43121337890625,
                            "Top": 0.876953125,
                        },
                        "Polygon": [
                            {"X": 0.43121337890625, "Y": 0.876953125},
                            {"X": 0.5679931640625, "Y": 0.876953125},
                            {"X": 0.5679931640625, "Y": 0.9072265625},
                            {"X": 0.43121337890625, "Y": 0.9072265625},
                        ],
                    },
                },
            },
            {
                "Timestamp": 1000,
                "TextDetection": {
                    "DetectedText": "Goodbye",
                    "Type": "WORD",
                    "Id": 1,
                    "ParentId": 0,
                    "Confidence": 99.7258529663086,
                    "Geometry": {
                        "BoundingBox": {
                            "Width": 0.0648193359375,
                            "Height": 0.0234375,
                            "Left": 0.43121337890625,
                            "Top": 0.876953125,
                        },
                        "Polygon": [
                            {"X": 0.43121337890625, "Y": 0.876953125},
                            {"X": 0.49603271484375, "Y": 0.876953125},
                            {"X": 0.49603271484375, "Y": 0.900390625},
                            {"X": 0.43121337890625, "Y": 0.900390625},
                        ],
                    },
                },
            },
            {
                "Timestamp": 1000,
                "TextDetection": {
                    "DetectedText": "world",
                    "Type": "WORD",
                    "Id": 2,
                    "ParentId": 0,
                    "Confidence": 98.22015380859375,
                    "Geometry": {
                        "BoundingBox": {
                            "Width": 0.0703125,
                            "Height": 0.0263671875,
                            "Left": 0.4976806640625,
                            "Top": 0.880859375,
                        },
                        "Polygon": [
                            {"X": 0.4976806640625, "Y": 0.880859375},
                            {"X": 0.5679931640625, "Y": 0.880859375},
                            {"X": 0.5679931640625, "Y": 0.9072265625},
                            {"X": 0.4976806640625, "Y": 0.9072265625},
                        ],
                    },
                },
            },
        ]

    def _face_matches(self) -> List[Dict[str, Any]]:
        return [
            {
                "Face": {
                    "BoundingBox": {
                        "Width": 0.5521978139877319,
                        "Top": 0.1203877404332161,
                        "Left": 0.23626373708248138,
                        "Height": 0.3126954436302185,
                    },
                    "Confidence": 99.98751068115234,
                    "Pose": {
                        "Yaw": -82.36799621582031,
                        "Roll": -62.13221740722656,
                        "Pitch": 0.8652129173278809,
                    },
                    "Quality": {
                        "Sharpness": 99.99880981445312,
                        "Brightness": 54.49755096435547,
                    },
                    "Landmarks": [
                        {
                            "Y": 0.2996366024017334,
                            "X": 0.41685718297958374,
                            "Type": "eyeLeft",
                        },
                        {
                            "Y": 0.2658946216106415,
                            "X": 0.4414493441581726,
                            "Type": "eyeRight",
                        },
                        {
                            "Y": 0.3465650677680969,
                            "X": 0.48636093735694885,
                            "Type": "nose",
                        },
                        {
                            "Y": 0.30935320258140564,
                            "X": 0.6251809000968933,
                            "Type": "mouthLeft",
                        },
                        {
                            "Y": 0.26942989230155945,
                            "X": 0.6454493403434753,
                            "Type": "mouthRight",
                        },
                    ],
                },
                "Similarity": 100.0,
            }
        ]

    def _unmatched_faces(self) -> List[Dict[str, Any]]:
        return [
            {
                "BoundingBox": {
                    "Width": 0.4890109896659851,
                    "Top": 0.6566604375839233,
                    "Left": 0.10989011079072952,
                    "Height": 0.278298944234848,
                },
                "Confidence": 99.99992370605469,
                "Pose": {
                    "Yaw": 51.51519012451172,
                    "Roll": -110.32493591308594,
                    "Pitch": -2.322134017944336,
                },
                "Quality": {
                    "Sharpness": 99.99671173095703,
                    "Brightness": 57.23163986206055,
                },
                "Landmarks": [
                    {
                        "Y": 0.8288310766220093,
                        "X": 0.3133862614631653,
                        "Type": "eyeLeft",
                    },
                    {
                        "Y": 0.7632885575294495,
                        "X": 0.28091415762901306,
                        "Type": "eyeRight",
                    },
                    {"Y": 0.7417283654212952, "X": 0.3631140887737274, "Type": "nose"},
                    {
                        "Y": 0.8081989884376526,
                        "X": 0.48565614223480225,
                        "Type": "mouthLeft",
                    },
                    {
                        "Y": 0.7548204660415649,
                        "X": 0.46090251207351685,
                        "Type": "mouthRight",
                    },
                ],
            }
        ]

    def source_image_face(self) -> Dict[str, Any]:
        return {
            "BoundingBox": {
                "Width": 0.5521978139877319,
                "Top": 0.1203877404332161,
                "Left": 0.23626373708248138,
                "Height": 0.3126954436302185,
            },
            "Confidence": 99.98751068115234,
        }

    def _mobile_phone_label(self) -> List[Dict[str, Any]]:
        return [
            {
                "Name": "Mobile Phone",
                "Parents": [{"Name": "Phone"}],
                "Aliases": [{"Name": "Cell Phone"}],
                "Categories": [{"Name": "Technology and Computing"}],
                "Confidence": 99.9364013671875,
                "Instances": [
                    {
                        "BoundingBox": {
                            "Width": 0.26779675483703613,
                            "Height": 0.8562285900115967,
                            "Left": 0.3604024350643158,
                            "Top": 0.09245597571134567,
                        },
                        "Confidence": 99.9364013671875,
                        "DominantColors": [
                            {
                                "Red": 120,
                                "Green": 137,
                                "Blue": 132,
                                "HexCode": "3A7432",
                                "SimplifiedColor": "red",
                                "CssColor": "fuscia",
                                "PixelPercentage": 40.10,
                            }
                        ],
                    }
                ],
            }
        ]

    def _image_properties(self) -> Dict[str, Any]:
        return {
            "ImageProperties": {
                "Quality": {
                    "Brightness": 40,
                    "Sharpness": 40,
                    "Contrast": 24,
                },
                "DominantColors": [
                    {
                        "Red": 120,
                        "Green": 137,
                        "Blue": 132,
                        "HexCode": "3A7432",
                        "SimplifiedColor": "red",
                        "CssColor": "fuscia",
                        "PixelPercentage": 40.10,
                    }
                ],
                "Foreground": {
                    "Quality": {
                        "Brightness": 40,
                        "Sharpness": 40,
                    },
                    "DominantColors": [
                        {
                            "Red": 200,
                            "Green": 137,
                            "Blue": 132,
                            "HexCode": "3A7432",
                            "CSSColor": "",
                            "SimplifiedColor": "red",
                            "PixelPercentage": 30.70,
                        }
                    ],
                },
                "Background": {
                    "Quality": {
                        "Brightness": 40,
                        "Sharpness": 40,
                    },
                    "DominantColors": [
                        {
                            "Red": 200,
                            "Green": 137,
                            "Blue": 132,
                            "HexCode": "3A7432",
                            "CSSColor": "",
                            "SimplifiedColor": "Red",
                            "PixelPercentage": 10.20,
                        }
                    ],
                },
            }
        }

    def _detect_text_text_detections(self) -> List[Dict[str, Any]]:
        return [
            {
                "Confidence": 99.35693359375,
                "DetectedText": "IT'S",
                "Geometry": {
                    "BoundingBox": {
                        "Height": 0.09988046437501907,
                        "Left": 0.6684935688972473,
                        "Top": 0.18226495385169983,
                        "Width": 0.1461552083492279,
                    },
                    "Polygon": [
                        {"X": 0.6684935688972473, "Y": 0.1838926374912262},
                        {"X": 0.8141663074493408, "Y": 0.18226495385169983},
                        {"X": 0.8146487474441528, "Y": 0.28051772713661194},
                        {"X": 0.6689760088920593, "Y": 0.2821454107761383},
                    ],
                },
                "Id": 0,
                "Type": "LINE",
            },
            {
                "Confidence": 99.6207275390625,
                "DetectedText": "MONDAY",
                "Geometry": {
                    "BoundingBox": {
                        "Height": 0.11442459374666214,
                        "Left": 0.5566731691360474,
                        "Top": 0.3525116443634033,
                        "Width": 0.39574965834617615,
                    },
                    "Polygon": [
                        {"X": 0.5566731691360474, "Y": 0.353712260723114},
                        {"X": 0.9522717595100403, "Y": 0.3525116443634033},
                        {"X": 0.9524227976799011, "Y": 0.4657355844974518},
                        {"X": 0.5568241477012634, "Y": 0.46693623065948486},
                    ],
                },
                "Id": 1,
                "Type": "LINE",
            },
            {
                "Confidence": 99.6160888671875,
                "DetectedText": "but keep",
                "Geometry": {
                    "BoundingBox": {
                        "Height": 0.08314694464206696,
                        "Left": 0.6398131847381592,
                        "Top": 0.5267938375473022,
                        "Width": 0.2021435648202896,
                    },
                    "Polygon": [
                        {"X": 0.640289306640625, "Y": 0.5267938375473022},
                        {"X": 0.8419567942619324, "Y": 0.5295097827911377},
                        {"X": 0.8414806723594666, "Y": 0.609940767288208},
                        {"X": 0.6398131847381592, "Y": 0.6072247624397278},
                    ],
                },
                "Id": 2,
                "Type": "LINE",
            },
            {
                "Confidence": 88.95134735107422,
                "DetectedText": "Smiling",
                "Geometry": {
                    "BoundingBox": {
                        "Height": 0.4326171875,
                        "Left": 0.46289217472076416,
                        "Top": 0.5634765625,
                        "Width": 0.5371078252792358,
                    },
                    "Polygon": [
                        {"X": 0.46289217472076416, "Y": 0.5634765625},
                        {"X": 1.0, "Y": 0.5634765625},
                        {"X": 1.0, "Y": 0.99609375},
                        {"X": 0.46289217472076416, "Y": 0.99609375},
                    ],
                },
                "Id": 3,
                "Type": "LINE",
            },
            {
                "Confidence": 99.35693359375,
                "DetectedText": "IT'S",
                "Geometry": {
                    "BoundingBox": {
                        "Height": 0.09988046437501907,
                        "Left": 0.6684935688972473,
                        "Top": 0.18226495385169983,
                        "Width": 0.1461552083492279,
                    },
                    "Polygon": [
                        {"X": 0.6684935688972473, "Y": 0.1838926374912262},
                        {"X": 0.8141663074493408, "Y": 0.18226495385169983},
                        {"X": 0.8146487474441528, "Y": 0.28051772713661194},
                        {"X": 0.6689760088920593, "Y": 0.2821454107761383},
                    ],
                },
                "Id": 4,
                "ParentId": 0,
                "Type": "WORD",
            },
            {
                "Confidence": 99.6207275390625,
                "DetectedText": "MONDAY",
                "Geometry": {
                    "BoundingBox": {
                        "Height": 0.11442466825246811,
                        "Left": 0.5566731691360474,
                        "Top": 0.35251158475875854,
                        "Width": 0.39574965834617615,
                    },
                    "Polygon": [
                        {"X": 0.5566731691360474, "Y": 0.3537122905254364},
                        {"X": 0.9522718787193298, "Y": 0.35251158475875854},
                        {"X": 0.9524227976799011, "Y": 0.4657355546951294},
                        {"X": 0.5568241477012634, "Y": 0.46693626046180725},
                    ],
                },
                "Id": 5,
                "ParentId": 1,
                "Type": "WORD",
            },
            {
                "Confidence": 99.96778869628906,
                "DetectedText": "but",
                "Geometry": {
                    "BoundingBox": {
                        "Height": 0.0625,
                        "Left": 0.6402802467346191,
                        "Top": 0.5283203125,
                        "Width": 0.08027780801057816,
                    },
                    "Polygon": [
                        {"X": 0.6402802467346191, "Y": 0.5283203125},
                        {"X": 0.7205580472946167, "Y": 0.5283203125},
                        {"X": 0.7205580472946167, "Y": 0.5908203125},
                        {"X": 0.6402802467346191, "Y": 0.5908203125},
                    ],
                },
                "Id": 6,
                "ParentId": 2,
                "Type": "WORD",
            },
            {
                "Confidence": 99.26438903808594,
                "DetectedText": "keep",
                "Geometry": {
                    "BoundingBox": {
                        "Height": 0.0818721204996109,
                        "Left": 0.7344760298728943,
                        "Top": 0.5280686020851135,
                        "Width": 0.10748066753149033,
                    },
                    "Polygon": [
                        {"X": 0.7349520921707153, "Y": 0.5280686020851135},
                        {"X": 0.8419566750526428, "Y": 0.5295097827911377},
                        {"X": 0.8414806127548218, "Y": 0.6099407076835632},
                        {"X": 0.7344760298728943, "Y": 0.6084995269775391},
                    ],
                },
                "Id": 7,
                "ParentId": 2,
                "Type": "WORD",
            },
            {
                "Confidence": 88.95134735107422,
                "DetectedText": "Smiling",
                "Geometry": {
                    "BoundingBox": {
                        "Height": 0.4326171875,
                        "Left": 0.46289217472076416,
                        "Top": 0.5634765625,
                        "Width": 0.5371078252792358,
                    },
                    "Polygon": [
                        {"X": 0.46289217472076416, "Y": 0.5634765625},
                        {"X": 1.0, "Y": 0.5634765625},
                        {"X": 1.0, "Y": 0.99609375},
                        {"X": 0.46289217472076416, "Y": 0.99609375},
                    ],
                },
                "Id": 8,
                "ParentId": 3,
                "Type": "WORD",
            },
        ]

    def _detect_custom_labels_detections(self) -> List[Dict[str, Any]]:
        return [
            {
                "Name": "MyLogo",
                "Confidence": 77.7729721069336,
                "Geometry": {
                    "BoundingBox": {
                        "Width": 0.198987677693367,
                        "Height": 0.31296101212501526,
                        "Left": 0.07924537360668182,
                        "Top": 0.4037395715713501,
                    }
                },
            }
        ]


rekognition_backends = BackendDict(RekognitionBackend, "rekognition")
