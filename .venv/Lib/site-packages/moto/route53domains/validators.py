import re
from datetime import datetime, timedelta, timezone
from enum import Enum
from ipaddress import IPv4Address, IPv6Address, ip_address
from typing import Any, Dict, List, Optional, Type

from moto.core.common_models import BaseModel
from moto.moto_api._internal import mock_random
from moto.route53domains.exceptions import UnsupportedTLDException

DOMAIN_OPERATION_STATUSES = (
    "SUBMITTED",
    "IN_PROGRESS",
    "ERROR",
    "SUCCESSFUL",
    "FAILED",
)

DOMAIN_OPERATION_TYPES = (
    "REGISTER_DOMAIN",
    "DELETE_DOMAIN",
    "TRANSFER_IN_DOMAIN",
    "UPDATE_DOMAIN_CONTACT",
    "UPDATE_NAMESERVER",
    "CHANGE_PRIVACY_PROTECTION",
    "DOMAIN_LOCK",
    "ENABLE_AUTORENEW",
    "DISABLE_AUTORENEW",
    "ADD_DNSSEC",
    "REMOVE_DNSSEC",
    "EXPIRE_DOMAIN",
    "TRANSFER_OUT_DOMAIN",
    "CHANGE_DOMAIN_OWNER",
    "RENEW_DOMAIN",
    "PUSH_DOMAIN",
    "INTERNAL_TRANSFER_OUT_DOMAIN",
    "INTERNAL_TRANSFER_IN_DOMAIN",
)

DOMAIN_OPERATION_STATUS_FLAGS = (
    "PENDING_ACCEPTANCE",
    "PENDING_CUSTOMER_ACTION",
    "PENDING_AUTHORIZATION",
    "PENDING_PAYMENT_VERIFICATION",
    "PENDING_SUPPORT_CASE",
)

DOMAIN_CONTACT_DETAIL_CONTACT_TYPES = (
    "PERSON",
    "COMPANY",
    "ASSOCIATION",
    "PUBLIC_BODY",
    "RESELLER",
    "ORGANIZATION",
)

DOMAIN_CONTACT_DETAIL_COUNTRY_CODES = (
    "AC",
    "AD",
    "AE",
    "AF",
    "AG",
    "AI",
    "AL",
    "AM",
    "AN",
    "AO",
    "AQ",
    "AR",
    "AS",
    "AT",
    "AU",
    "AW",
    "AX",
    "AZ",
    "BA",
    "BB",
    "BD",
    "BE",
    "BF",
    "BG",
    "BH",
    "BI",
    "BJ",
    "BL",
    "BM",
    "BN",
    "BO",
    "BQ",
    "BR",
    "BS",
    "BT",
    "BV",
    "BW",
    "BY",
    "BZ",
    "CA",
    "CC",
    "CD",
    "CF",
    "CG",
    "CH",
    "CI",
    "CK",
    "CL",
    "CM",
    "CN",
    "CO",
    "CR",
    "CU",
    "CV",
    "CW",
    "CX",
    "CY",
    "CZ",
    "DE",
    "DJ",
    "DK",
    "DM",
    "DO",
    "DZ",
    "EC",
    "EE",
    "EG",
    "EH",
    "ER",
    "ES",
    "ET",
    "FI",
    "FJ",
    "FK",
    "FM",
    "FO",
    "FR",
    "GA",
    "GB",
    "GD",
    "GE",
    "GF",
    "GG",
    "GH",
    "GI",
    "GL",
    "GM",
    "GN",
    "GP",
    "GQ",
    "GR",
    "GS",
    "GT",
    "GU",
    "GW",
    "GY",
    "HK",
    "HM",
    "HN",
    "HR",
    "HT",
    "HU",
    "ID",
    "IE",
    "IL",
    "IM",
    "IN",
    "IO",
    "IQ",
    "IR",
    "IS",
    "IT",
    "JE",
    "JM",
    "JO",
    "JP",
    "KE",
    "KG",
    "KH",
    "KI",
    "KM",
    "KN",
    "KP",
    "KR",
    "KW",
    "KY",
    "KZ",
    "LA",
    "LB",
    "LC",
    "LI",
    "LK",
    "LR",
    "LS",
    "LT",
    "LU",
    "LV",
    "LY",
    "MA",
    "MC",
    "MD",
    "ME",
    "MF",
    "MG",
    "MH",
    "MK",
    "ML",
    "MM",
    "MN",
    "MO",
    "MP",
    "MQ",
    "MR",
    "MS",
    "MT",
    "MU",
    "MV",
    "MW",
    "MX",
    "MY",
    "MZ",
    "NA",
    "NC",
    "NE",
    "NF",
    "NG",
    "NI",
    "NL",
    "NO",
    "NP",
    "NR",
    "NU",
    "NZ",
    "OM",
    "PA",
    "PE",
    "PF",
    "PG",
    "PH",
    "PK",
    "PL",
    "PM",
    "PN",
    "PR",
    "PS",
    "PT",
    "PW",
    "PY",
    "QA",
    "RE",
    "RO",
    "RS",
    "RU",
    "RW",
    "SA",
    "SB",
    "SC",
    "SD",
    "SE",
    "SG",
    "SH",
    "SI",
    "SJ",
    "SK",
    "SL",
    "SM",
    "SN",
    "SO",
    "SR",
    "SS",
    "ST",
    "SV",
    "SX",
    "SY",
    "SZ",
    "TC",
    "TD",
    "TF",
    "TG",
    "TH",
    "TJ",
    "TK",
    "TL",
    "TM",
    "TN",
    "TO",
    "TP",
    "TR",
    "TT",
    "TV",
    "TW",
    "TZ",
    "UA",
    "UG",
    "US",
    "UY",
    "UZ",
    "VA",
    "VC",
    "VE",
    "VG",
    "VI",
    "VN",
    "VU",
    "WF",
    "WS",
    "YE",
    "YT",
    "ZA",
    "ZM",
    "ZW",
)

# List of supported top-level domains that you can register with Amazon Route53
AWS_SUPPORTED_TLDS = (
    "ac",
    "academy",
    "accountants",
    "actor",
    "adult",
    "agency",
    "airforce",
    "apartments",
    "associates",
    "auction",
    "audio",
    "band",
    "bargains",
    "bet",
    "bike",
    "bingo",
    "biz",
    "black",
    "blue",
    "boutique",
    "builders",
    "business",
    "buzz",
    "cab",
    "cafe",
    "camera",
    "camp",
    "capital",
    "cards",
    "care",
    "careers",
    "cash",
    "casino",
    "catering",
    "cc",
    "center",
    "ceo",
    "chat",
    "cheap",
    "church",
    "city",
    "claims",
    "cleaning",
    "click",
    "clinic",
    "clothing",
    "cloud",
    "club",
    "coach",
    "codes",
    "coffee",
    "college",
    "com",
    "community",
    "company",
    "computer",
    "condos",
    "construction",
    "consulting",
    "contractors",
    "cool",
    "coupons",
    "credit",
    "creditcard",
    "cruises",
    "dance",
    "dating",
    "deals",
    "degree",
    "delivery",
    "democrat",
    "dental",
    "diamonds",
    "diet",
    "digital",
    "direct",
    "directory",
    "discount",
    "dog",
    "domains",
    "education",
    "email",
    "energy",
    "engineering",
    "enterprises",
    "equipment",
    "estate",
    "events",
    "exchange",
    "expert",
    "exposed",
    "express",
    "fail",
    "farm",
    "finance",
    "financial",
    "fish",
    "fitness",
    "flights",
    "florist",
    "flowers",
    "fm",
    "football",
    "forsale",
    "foundation",
    "fund",
    "furniture",
    "futbol",
    "fyi",
    "gallery",
    "games",
    "gift",
    "gifts",
    "gives",
    "glass",
    "global",
    "gmbh",
    "gold",
    "golf",
    "graphics",
    "gratis",
    "green",
    "gripe",
    "group",
    "guide",
    "guitars",
    "guru",
    "haus",
    "healthcare",
    "help",
    "hiv",
    "hockey",
    "holdings",
    "holiday",
    "host",
    "hosting",
    "house",
    "im",
    "immo",
    "immobilien",
    "industries",
    "info",
    "ink",
    "institute",
    "insure",
    "international",
    "investments",
    "io",
    "irish",
    "jewelry",
    "juegos",
    "kaufen",
    "kim",
    "kitchen",
    "kiwi",
    "land",
    "lease",
    "legal",
    "lgbt",
    "life",
    "lighting",
    "limited",
    "limo",
    "link",
    "live",
    "loan",
    "loans",
    "lol",
    "maison",
    "management",
    "marketing",
    "mba",
    "media",
    "memorial",
    "mobi",
    "moda",
    "money",
    "mortgage",
    "movie",
    "name",
    "net",
    "network",
    "news",
    "ninja",
    "onl",
    "online",
    "org",
    "partners",
    "parts",
    "photo",
    "photography",
    "photos",
    "pics",
    "pictures",
    "pink",
    "pizza",
    "place",
    "plumbing",
    "plus",
    "poker",
    "porn",
    "press",
    "pro",
    "productions",
    "properties",
    "property",
    "pub",
    "qpon",
    "recipes",
    "red",
    "reise",
    "reisen",
    "rentals",
    "repair",
    "report",
    "republican",
    "restaurant",
    "reviews",
    "rip",
    "rocks",
    "run",
    "sale",
    "sarl",
    "school",
    "schule",
    "services",
    "sex",
    "sexy",
    "shiksha",
    "shoes",
    "show",
    "singles",
    "site",
    "soccer",
    "social",
    "solar",
    "solutions",
    "space",
    "store",
    "studio",
    "style",
    "sucks",
    "supplies",
    "supply",
    "support",
    "surgery",
    "systems",
    "tattoo",
    "tax",
    "taxi",
    "team",
    "tech",
    "technology",
    "tennis",
    "theater",
    "tienda",
    "tips",
    "tires",
    "today",
    "tools",
    "tours",
    "town",
    "toys",
    "trade",
    "training",
    "tv",
    "university",
    "uno",
    "vacations",
    "vegas",
    "ventures",
    "vg",
    "viajes",
    "video",
    "villas",
    "vision",
    "voyage",
    "watch",
    "website",
    "wedding",
    "wiki",
    "wine",
    "works",
    "world",
    "wtf",
    "xyz",
    "zone",
)

VALID_DOMAIN_REGEX = re.compile(
    r"(?:[a-z0-9](?:[a-z0-9-]{0,61}[a-z0-9])?\.)+[a-z0-9][a-z0-9-]{0,61}[a-z0-9]"
)
PHONE_NUMBER_REGEX = re.compile(r"\+\d*\.\d+$")


class DomainFilterField(str, Enum):
    DOMAIN_NAME = "DomainName"
    EXPIRY = "Expiry"


class DomainSortOrder(str, Enum):
    ASCENDING = "ASC"
    DESCENDING = "DES"


class DomainFilterOperator(str, Enum):
    LE = "LE"
    GE = "GE"
    BEGINS_WITH = "BEGINS_WITH"


def is_valid_enum(value: Any, enum_cls: Type[Enum]) -> bool:
    try:
        enum_cls(value)
        return True
    except ValueError:
        return False


class ValidationException(Exception):
    def __init__(self, errors: List[str]):
        super().__init__("\n\t".join(errors))
        self.errors = errors


class Route53DomainsOperation(BaseModel):
    def __init__(
        self,
        id_: str,
        domain_name: str,
        status: str,
        type_: str,
        submitted_date: datetime,
        last_updated_date: datetime,
        message: Optional[str] = None,
        status_flag: Optional[str] = None,
    ):
        self.id = id_
        self.domain_name = domain_name
        self.status = status
        self.type = type_
        self.submitted_date = submitted_date
        self.last_updated_date = last_updated_date
        self.message = message
        self.status_flag = status_flag

    @classmethod
    def validate(  # type: ignore[misc,no-untyped-def]
        cls,
        domain_name: str,
        status: str,
        type_: str,
        message: Optional[str] = None,
        status_flag: Optional[str] = None,
    ):
        id_ = str(mock_random.uuid4())
        submitted_date = datetime.now(timezone.utc)
        last_updated_date = datetime.now(timezone.utc)

        return cls(
            id_,
            domain_name,
            status,
            type_,
            submitted_date,
            last_updated_date,
            message,
            status_flag,
        )

    def to_json(self) -> Dict[str, Any]:
        d = {
            "OperationId": self.id,
            "Status": self.status,
            "StatusFlag": self.status_flag,
            "DomainName": self.domain_name,
            "LastUpdatedDate": self.last_updated_date.timestamp(),
            "SubmittedDate": self.submitted_date.timestamp(),
            "Type": self.type,
        }

        if self.message:
            d["Message"] = self.message

        return d


class Route53DomainsContactDetail(BaseModel):
    def __init__(
        self,
        address_line_1: Optional[str] = None,
        address_line_2: Optional[str] = None,
        city: Optional[str] = None,
        contact_type: Optional[str] = None,
        country_code: Optional[str] = None,
        email: Optional[str] = None,
        extra_params: Optional[List[Dict[str, Any]]] = None,
        fax: Optional[str] = None,
        first_name: Optional[str] = None,
        last_name: Optional[str] = None,
        organization_name: Optional[str] = None,
        phone_number: Optional[str] = None,
        state: Optional[str] = None,
        zip_code: Optional[str] = None,
    ):
        super().__init__()
        self.address_line_1 = address_line_1
        self.address_line_2 = address_line_2
        self.city = city
        self.contact_type = contact_type
        self.country_code = country_code
        self.email = email
        self.extra_params = extra_params
        self.fax = fax
        self.first_name = first_name
        self.last_name = last_name
        self.organization_name = organization_name
        self.phone_number = phone_number
        self.state = state
        self.zip_code = zip_code

    @classmethod
    def validate(  # type: ignore[misc, no-untyped-def]
        cls,
        address_line_1: Optional[str] = None,
        address_line_2: Optional[str] = None,
        city: Optional[str] = None,
        contact_type: Optional[str] = None,
        country_code: Optional[str] = None,
        email: Optional[str] = None,
        extra_params: Optional[List[Dict[str, Any]]] = None,
        fax: Optional[str] = None,
        first_name: Optional[str] = None,
        last_name: Optional[str] = None,
        organization_name: Optional[str] = None,
        phone_number: Optional[str] = None,
        state: Optional[str] = None,
        zip_code: Optional[str] = None,
    ):
        input_errors: List[str] = []

        cls.__validate_str_len(address_line_1, "AddressLine1", 255, input_errors)
        cls.__validate_str_len(address_line_2, "AddressLine2", 255, input_errors)
        cls.__validate_str_len(city, "City", 255, input_errors)
        cls.__validate_str_len(email, "Email", 255, input_errors)
        cls.__validate_str_len(fax, "Fax", 255, input_errors)
        cls.__validate_str_len(first_name, "FirstName", 255, input_errors)
        cls.__validate_str_len(last_name, "LastName", 255, input_errors)
        cls.__validate_str_len(state, "State", 255, input_errors)
        cls.__validate_str_len(zip_code, "ZipCode", 255, input_errors)

        if contact_type:
            if contact_type not in DOMAIN_CONTACT_DETAIL_CONTACT_TYPES:
                input_errors.append(f"Invalid contact type {contact_type}")
            else:
                if contact_type != "PERSON" and not organization_name:
                    input_errors.append(
                        "Must supply OrganizationName when ContactType is not PERSON"
                    )

        if country_code and country_code not in DOMAIN_CONTACT_DETAIL_COUNTRY_CODES:
            input_errors.append(f"CountryCode {country_code} is invalid")

        if phone_number and not PHONE_NUMBER_REGEX.match(phone_number):
            input_errors.append("PhoneNumber is in an invalid format")

        if input_errors:
            raise ValidationException(input_errors)

        return cls(
            address_line_1,
            address_line_2,
            city,
            contact_type,
            country_code,
            email,
            extra_params,
            fax,
            first_name,
            last_name,
            organization_name,
            phone_number,
            state,
            zip_code,
        )

    @classmethod
    def validate_dict(cls, d: Dict[str, Any]):  # type: ignore[misc, no-untyped-def]
        address_line_1 = d.get("AddressLine1")
        address_line_2 = d.get("AddressLine2")
        city = d.get("City")
        contact_type = d.get("ContactType")
        country_code = d.get("CountryCode")
        email = d.get("Email")
        extra_params = d.get("ExtraParams")
        fax = d.get("Fax")
        first_name = d.get("FirstName")
        last_name = d.get("LastName")
        organization_name = d.get("OrganizationName")
        phone_number = d.get("PhoneNumber")
        state = d.get("State")
        zip_code = d.get("ZipCode")
        return cls.validate(
            address_line_1=address_line_1,
            address_line_2=address_line_2,
            city=city,
            contact_type=contact_type,
            country_code=country_code,
            email=email,
            extra_params=extra_params,
            fax=fax,
            first_name=first_name,
            last_name=last_name,
            organization_name=organization_name,
            phone_number=phone_number,
            state=state,
            zip_code=zip_code,
        )

    @staticmethod
    def __validate_str_len(
        value: Optional[str], field_name: str, max_len: int, input_errors: List[str]
    ) -> None:
        if value and len(value) > max_len:
            input_errors.append(f"Length of {field_name} is more than {max_len}")

    def to_json(self) -> Dict[str, Any]:
        d = {
            "FirstName": self.first_name,
            "LastName": self.last_name,
            "ContactType": self.contact_type,
            "OrganizationName": self.organization_name,
            "AddressLine1": self.address_line_1,
            "AddressLine2": self.address_line_2,
            "City": self.city,
            "State": self.state,
            "CountryCode": self.country_code,
            "ZipCode": self.zip_code,
            "PhoneNumber": self.phone_number,
            "Email": self.email,
            "Fax": self.fax,
            "ExtraParams": self.extra_params,
        }

        return {key: value for key, value in d.items() if value is not None}


class NameServer:
    def __init__(self, name: str, glue_ips: List[str]):
        self.name = name
        self.glue_ips = glue_ips

    @classmethod
    def validate(cls, name: str, glue_ips: Optional[List[str]] = None):  # type: ignore[misc,no-untyped-def]
        glue_ips = glue_ips or []
        input_errors: List[str] = []

        if not VALID_DOMAIN_REGEX.match(name):
            input_errors.append(f"{name} is not a valid host name")

        num_ipv4_addresses = 0
        num_ipv6_addresses = 0
        for ip in glue_ips:
            try:
                address = ip_address(ip)
                if isinstance(address, IPv4Address):
                    num_ipv4_addresses += 1
                elif isinstance(address, IPv6Address):
                    num_ipv6_addresses += 1
            except ValueError:
                input_errors.append(f"{ip} is not a valid IP address")

        if num_ipv4_addresses > 1:
            input_errors.append("GlueIps list must include only 1 IPv4 address")
        if num_ipv6_addresses > 1:
            input_errors.append("GlueIps list must include only 1 IPv6 address")

        if input_errors:
            raise ValidationException(input_errors)

        return cls(name, glue_ips)

    @classmethod
    def validate_dict(cls, data: Dict[str, Any]):  # type: ignore[misc,no-untyped-def]
        name = data.get("Name")
        glue_ips = data.get("GlueIps")
        return cls.validate(name, glue_ips)  # type: ignore[arg-type]

    def to_json(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {"Name": self.name}
        if self.glue_ips:
            d["GlueIps"] = self.glue_ips
        return d


class Route53Domain(BaseModel):
    def __init__(
        self,
        domain_name: str,
        nameservers: List[NameServer],
        auto_renew: bool,
        admin_contact: Route53DomainsContactDetail,
        registrant_contact: Route53DomainsContactDetail,
        tech_contact: Route53DomainsContactDetail,
        admin_privacy: bool,
        registrant_privacy: bool,
        tech_privacy: bool,
        registrar_name: str,
        whois_server: str,
        registrar_url: str,
        abuse_contact_email: str,
        abuse_contact_phone: str,
        registry_domain_id: str,
        creation_date: datetime,
        updated_date: datetime,
        expiration_date: datetime,
        reseller: str,
        status_list: List[str],
        dns_sec_keys: List[Dict[str, Any]],
        extra_params: List[Dict[str, Any]],
    ):
        self.domain_name = domain_name
        self.nameservers = nameservers
        self.auto_renew = auto_renew
        self.admin_contact = admin_contact
        self.registrant_contact = registrant_contact
        self.tech_contact = tech_contact
        self.admin_privacy = admin_privacy
        self.registrant_privacy = registrant_privacy
        self.tech_privacy = tech_privacy
        self.registrar_name = registrar_name
        self.whois_server = whois_server
        self.registrar_url = registrar_url
        self.abuse_contact_email = abuse_contact_email
        self.abuse_contact_phone = abuse_contact_phone
        self.registry_domain_id = registry_domain_id
        self.creation_date = creation_date
        self.updated_date = updated_date
        self.expiration_date = expiration_date
        self.reseller = reseller
        self.status_list = status_list
        self.dns_sec_keys = dns_sec_keys
        self.extra_params = extra_params

    @classmethod
    def validate(  # type: ignore[misc,no-untyped-def]
        cls,
        domain_name: str,
        admin_contact: Route53DomainsContactDetail,
        registrant_contact: Route53DomainsContactDetail,
        tech_contact: Route53DomainsContactDetail,
        nameservers: Optional[List[Dict[str, Any]]] = None,
        auto_renew: bool = True,
        admin_privacy: bool = True,
        registrant_privacy: bool = True,
        tech_privacy: bool = True,
        registrar_name: Optional[str] = None,
        whois_server: Optional[str] = None,
        registrar_url: Optional[str] = None,
        abuse_contact_email: Optional[str] = None,
        abuse_contact_phone: Optional[str] = None,
        registry_domain_id: Optional[str] = None,
        expiration_date: Optional[datetime] = None,
        reseller: Optional[str] = None,
        dns_sec_keys: Optional[List[Dict[str, Any]]] = None,
        extra_params: Optional[List[Dict[str, Any]]] = None,
    ):
        input_errors: List[str] = []

        cls.validate_domain_name(domain_name, input_errors)

        nameservers = nameservers or []
        try:
            nameservers = [
                NameServer.validate_dict(nameserver) for nameserver in nameservers
            ] or [
                NameServer.validate(name="ns-2048.awscdn-64.net"),
                NameServer.validate(name="ns-2051.awscdn-67.net"),
                NameServer.validate(name="ns-2050.awscdn-66.net"),
                NameServer.validate(name="ns-2049.awscdn-65.net"),
            ]
        except ValidationException as e:
            input_errors += e.errors

        creation_date = datetime.now(timezone.utc)
        updated_date = datetime.now(timezone.utc)
        expiration_date = expiration_date or datetime.now(timezone.utc) + timedelta(
            days=365 * 10
        )
        registrar_name = registrar_name or "GANDI SAS"
        whois_server = whois_server or "whois.gandi.net"
        registrar_url = registrar_url or "http://www.gandi.net"
        abuse_contact_email = abuse_contact_email or "abuse@support.gandi.net"
        status_list = ["SUCCEEDED"]

        time_until_expiration = expiration_date - datetime.now(timezone.utc)
        if time_until_expiration < timedelta(
            days=365
        ) or time_until_expiration > timedelta(days=365 * 10):
            input_errors.append(
                "ExpirationDate must by between 1 and 10 years from now"
            )

        if input_errors:
            raise ValidationException(input_errors)

        return cls(
            domain_name=domain_name,
            nameservers=nameservers,  # type: ignore[arg-type]
            auto_renew=auto_renew,
            admin_contact=admin_contact,
            registrant_contact=registrant_contact,
            tech_contact=tech_contact,
            admin_privacy=admin_privacy,
            registrant_privacy=registrant_privacy,
            tech_privacy=tech_privacy,
            registrar_name=registrar_name,
            whois_server=whois_server,
            registrar_url=registrar_url,
            abuse_contact_email=abuse_contact_email,
            abuse_contact_phone=abuse_contact_phone,  # type: ignore[arg-type]
            registry_domain_id=registry_domain_id,  # type: ignore[arg-type]
            creation_date=creation_date,
            updated_date=updated_date,
            expiration_date=expiration_date,
            reseller=reseller,  # type: ignore[arg-type]
            status_list=status_list,
            dns_sec_keys=dns_sec_keys,  # type: ignore[arg-type]
            extra_params=extra_params,  # type: ignore[arg-type]
        )

    @staticmethod
    def validate_domain_name(domain_name: str, input_errors: List[str]) -> None:
        if not VALID_DOMAIN_REGEX.match(domain_name):
            input_errors.append("Invalid domain name")
            return

        tld = domain_name.split(".")[-1]
        if tld not in AWS_SUPPORTED_TLDS:
            raise UnsupportedTLDException(tld)

    def to_json(self) -> Dict[str, Any]:
        return {
            "DomainName": self.domain_name,
            "Nameservers": [nameserver.to_json() for nameserver in self.nameservers],
            "AutoRenew": self.auto_renew,
            "AdminContact": self.admin_contact.to_json(),
            "RegistrantContact": self.registrant_contact.to_json(),
            "TechContact": self.tech_contact.to_json(),
            "AdminPrivacy": self.admin_privacy,
            "RegistrantPrivacy": self.registrant_privacy,
            "TechPrivacy": self.tech_privacy,
            "RegistrarName": self.registrar_name,
            "WhoIsServer": self.whois_server,
            "RegistrarUrl": self.registrar_url,
            "AbuseContactEmail": self.abuse_contact_email,
            "AbuseContactPhone": self.abuse_contact_phone,
            "RegistryDomainId": "",
            "CreationDate": self.creation_date.timestamp(),
            "UpdateDate": self.updated_date.timestamp(),
            "ExpirationDate": self.expiration_date.timestamp(),
            "Reseller": self.reseller,
            "DnsSec": "",
            "StatusList": self.status_list,
            "DnsSecKeys": self.dns_sec_keys,
            "BillingContact": self.admin_contact.to_json(),
        }


class DomainsFilter:
    def __init__(
        self, name: DomainFilterField, operator: DomainFilterOperator, values: List[str]
    ):
        self.name: DomainFilterField = name
        self.operator: DomainFilterOperator = operator
        self.values = values

    def filter(self, domain: Route53Domain) -> bool:
        if self.name == DomainFilterField.DOMAIN_NAME:
            return self.__filter_by_domain_name(domain)
        return self.__filter_by_expiry_date(domain)

    def __filter_by_domain_name(self, domain: Route53Domain) -> bool:
        return any([domain.domain_name.startswith(value) for value in self.values])

    def __filter_by_expiry_date(self, domain: Route53Domain) -> bool:
        return (
            any(
                [
                    value
                    for value in self.values
                    if domain.expiration_date
                    >= datetime.fromtimestamp(float(value), tz=timezone.utc)
                ]
            )
            if self.operator == DomainFilterOperator.GE
            else any(
                [
                    value
                    for value in self.values
                    if domain.expiration_date
                    <= datetime.fromtimestamp(float(value), tz=timezone.utc)
                ]
            )
        )

    @classmethod
    def validate(cls, name: str, operator: str, values: List[str]):  # type: ignore[misc, no-untyped-def]
        input_errors: List[str] = []

        if not is_valid_enum(name, DomainFilterField):
            input_errors.append(f"Cannot filter by field {name}")

        if not is_valid_enum(operator, DomainFilterOperator):
            input_errors.append(f"Invalid filter operator {operator}")

        if len(values) != 1:
            input_errors.append("Multiple filter values are not currently supported")

        if (
            name == DomainFilterField.DOMAIN_NAME
            and operator != DomainFilterOperator.BEGINS_WITH
        ):
            input_errors.append(
                f"Operator {operator} cannot be used with the DomainName filter"
            )

        if name == DomainFilterField.EXPIRY and operator not in (
            DomainFilterOperator.GE,
            DomainFilterOperator.LE,
        ):
            input_errors.append(
                f"Operator {operator} cannot be used with the Expiry filter"
            )

        if input_errors:
            raise ValidationException(input_errors)

        return cls(
            name=DomainFilterField(name),
            operator=DomainFilterOperator(operator),
            values=values,
        )

    @classmethod
    def validate_dict(cls, data: Dict[str, Any]):  # type: ignore[misc,no-untyped-def]
        name = data.get("Name")
        operator = data.get("Operator")
        values = data.get("Values")
        return cls.validate(name=name, operator=operator, values=values)  # type: ignore[arg-type]


class DomainsSortCondition:
    def __init__(self, name: DomainFilterField, sort_order: DomainSortOrder):
        self.name: DomainFilterField = name
        self.sort_order: DomainSortOrder = sort_order

    @classmethod
    def validate(cls, name: str, sort_order: str):  # type: ignore[misc,no-untyped-def]
        input_errors: List[str] = []
        if not is_valid_enum(name, DomainFilterField):
            input_errors.append(f"Cannot sort by field {name}")

        if not is_valid_enum(sort_order, DomainSortOrder):
            input_errors.append(f"Invalid sort order {sort_order}")

        if input_errors:
            raise ValidationException(input_errors)

        return cls(name=DomainFilterField(name), sort_order=DomainSortOrder(sort_order))

    @classmethod
    def validate_dict(cls, data: Dict[str, Any]):  # type: ignore[misc,no-untyped-def]
        name = data.get("Name")
        sort_order = data.get("SortOrder")
        return cls.validate(name=name, sort_order=sort_order)  # type: ignore[arg-type]
