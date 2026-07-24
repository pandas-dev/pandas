from typing import Final

class CreditCardNumbers:
    class CardTypeIndicators:
        Business: Final = "4229989800000003"
        Commercial: Final = "4111111111131010"
        Consumer: Final = "4229989700000004"
        Corporate: Final = "4229989100000000"
        DurbinRegulated: Final = "4111161010101010"
        Debit: Final = "4117101010101010"
        Healthcare: Final = "4111111510101010"
        Payroll: Final = "4111111114101010"
        Prepaid: Final = "4111111111111210"
        PrepaidReloadable: Final = "4229989900000002"
        Purchase: Final = "4229989500000006"
        IssuingBank: Final = "4111111141010101"
        CountryOfIssuance: Final = "4111111111121102"
        No: Final = "4111111111310101"
        Unknown: Final = "4111111111112101"

    Maestro: Final = "6304000000000000"
    MasterCard: Final = "5555555555554444"
    MasterCardInternational: Final = "5105105105105100"
    Visa: Final = "4012888888881881"
    VisaInternational: Final = "4009348888881881"
    VisaPrepaid: Final = "4500600000000061"
    Discover: Final = "6011111111111117"
    Elo: Final = "5066991111111118"
    Hiper: Final = "6370950000000005"
    Hipercard: Final = "6062820524845321"
    Amex: Final = "378734493671000"

    class FailsSandboxVerification:
        AmEx: Final = "378734493671000"
        Discover: Final = "6011000990139424"
        MasterCard: Final = "5105105105105100"
        Visa: Final = "4000111111111115"

    class AmexPayWithPoints:
        Success: Final = "371260714673002"
        IneligibleCard: Final = "378267515471109"
        InsufficientPoints: Final = "371544868764018"

    class Disputes:
        Chargeback: Final = "4023898493988028"
