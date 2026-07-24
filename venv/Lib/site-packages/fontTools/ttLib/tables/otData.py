from .otDataSchema import FieldSpec

otData = [
    #
    # common
    #
    ("LookupOrder", []),
    (
        "ScriptList",
        [
            FieldSpec("uint16", "ScriptCount", description="Number of ScriptRecords"),
            FieldSpec(
                "struct",
                "ScriptRecord",
                repeat="ScriptCount",
                aux=0,
                description="Array of ScriptRecords -listed alphabetically by ScriptTag",
            ),
        ],
    ),
    (
        "ScriptRecord",
        [
            FieldSpec("Tag", "ScriptTag", description="4-byte ScriptTag identifier"),
            FieldSpec(
                "Offset",
                "Script",
                description="Offset to Script table-from beginning of ScriptList",
            ),
        ],
    ),
    (
        "Script",
        [
            FieldSpec(
                "Offset",
                "DefaultLangSys",
                description="Offset to DefaultLangSys table-from beginning of Script table-may be NULL",
            ),
            FieldSpec(
                "uint16",
                "LangSysCount",
                description="Number of LangSysRecords for this script-excluding the DefaultLangSys",
            ),
            FieldSpec(
                "struct",
                "LangSysRecord",
                repeat="LangSysCount",
                aux=0,
                description="Array of LangSysRecords-listed alphabetically by LangSysTag",
            ),
        ],
    ),
    (
        "LangSysRecord",
        [
            FieldSpec("Tag", "LangSysTag", description="4-byte LangSysTag identifier"),
            FieldSpec(
                "Offset",
                "LangSys",
                description="Offset to LangSys table-from beginning of Script table",
            ),
        ],
    ),
    (
        "LangSys",
        [
            FieldSpec(
                "Offset",
                "LookupOrder",
                description="= NULL (reserved for an offset to a reordering table)",
            ),
            FieldSpec(
                "uint16",
                "ReqFeatureIndex",
                description="Index of a feature required for this language system- if no required features = 0xFFFF",
            ),
            FieldSpec(
                "uint16",
                "FeatureCount",
                description="Number of FeatureIndex values for this language system-excludes the required feature",
            ),
            FieldSpec(
                "uint16",
                "FeatureIndex",
                repeat="FeatureCount",
                aux=0,
                description="Array of indices into the FeatureList-in arbitrary order",
            ),
        ],
    ),
    (
        "FeatureList",
        [
            FieldSpec(
                "uint16",
                "FeatureCount",
                description="Number of FeatureRecords in this table",
            ),
            FieldSpec(
                "struct",
                "FeatureRecord",
                repeat="FeatureCount",
                aux=0,
                description="Array of FeatureRecords-zero-based (first feature has FeatureIndex = 0)-listed alphabetically by FeatureTag",
            ),
        ],
    ),
    (
        "FeatureRecord",
        [
            FieldSpec(
                "Tag", "FeatureTag", description="4-byte feature identification tag"
            ),
            FieldSpec(
                "Offset",
                "Feature",
                description="Offset to Feature table-from beginning of FeatureList",
            ),
        ],
    ),
    (
        "Feature",
        [
            FieldSpec(
                "Offset",
                "FeatureParams",
                description="= NULL (reserved for offset to FeatureParams)",
            ),
            FieldSpec(
                "uint16",
                "LookupCount",
                description="Number of LookupList indices for this feature",
            ),
            FieldSpec(
                "uint16",
                "LookupListIndex",
                repeat="LookupCount",
                aux=0,
                description="Array of LookupList indices for this feature -zero-based (first lookup is LookupListIndex = 0)",
            ),
        ],
    ),
    ("FeatureParams", []),
    (
        "FeatureParamsSize",
        [
            FieldSpec(
                "DeciPoints",
                "DesignSize",
                description="The design size in 720/inch units (decipoints).",
            ),
            FieldSpec(
                "uint16",
                "SubfamilyID",
                description="Serves as an identifier that associates fonts in a subfamily.",
            ),
            FieldSpec("NameID", "SubfamilyNameID", description="Subfamily NameID."),
            FieldSpec(
                "DeciPoints",
                "RangeStart",
                description="Small end of recommended usage range (exclusive) in 720/inch units.",
            ),
            FieldSpec(
                "DeciPoints",
                "RangeEnd",
                description="Large end of recommended usage range (inclusive) in 720/inch units.",
            ),
        ],
    ),
    (
        "FeatureParamsStylisticSet",
        [
            FieldSpec("uint16", "Version", description="Set to 0."),
            FieldSpec("NameID", "UINameID", description="UI NameID."),
        ],
    ),
    (
        "FeatureParamsCharacterVariants",
        [
            FieldSpec("uint16", "Format", description="Set to 0."),
            FieldSpec(
                "NameID", "FeatUILabelNameID", description="Feature UI label NameID."
            ),
            FieldSpec(
                "NameID",
                "FeatUITooltipTextNameID",
                description="Feature UI tooltip text NameID.",
            ),
            FieldSpec("NameID", "SampleTextNameID", description="Sample text NameID."),
            FieldSpec(
                "uint16",
                "NumNamedParameters",
                description="Number of named parameters.",
            ),
            FieldSpec(
                "NameID",
                "FirstParamUILabelNameID",
                description="First NameID of UI feature parameters.",
            ),
            FieldSpec(
                "uint16",
                "CharCount",
                description="Count of characters this feature provides glyph variants for.",
            ),
            FieldSpec(
                "uint24",
                "Character",
                repeat="CharCount",
                aux=0,
                description="Unicode characters for which this feature provides glyph variants.",
            ),
        ],
    ),
    (
        "LookupList",
        [
            FieldSpec(
                "uint16", "LookupCount", description="Number of lookups in this table"
            ),
            FieldSpec(
                "Offset",
                "Lookup",
                repeat="LookupCount",
                aux=0,
                description="Array of offsets to Lookup tables-from beginning of LookupList -zero based (first lookup is Lookup index = 0)",
            ),
        ],
    ),
    (
        "Lookup",
        [
            FieldSpec(
                "uint16",
                "LookupType",
                description="Different enumerations for GSUB and GPOS",
            ),
            FieldSpec("LookupFlag", "LookupFlag", description="Lookup qualifiers"),
            FieldSpec(
                "uint16",
                "SubTableCount",
                description="Number of SubTables for this lookup",
            ),
            FieldSpec(
                "Offset",
                "SubTable",
                repeat="SubTableCount",
                aux=0,
                description="Array of offsets to SubTables-from beginning of Lookup table",
            ),
            FieldSpec(
                "uint16",
                "MarkFilteringSet",
                aux="LookupFlag & 0x0010",
                description="If set, indicates that the lookup table structure is followed by a MarkFilteringSet field. The layout engine skips over all mark glyphs not in the mark filtering set indicated.",
            ),
        ],
    ),
    (
        "CoverageFormat1",
        [
            FieldSpec(
                "uint16", "CoverageFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "uint16", "GlyphCount", description="Number of glyphs in the GlyphArray"
            ),
            FieldSpec(
                "GlyphID",
                "GlyphArray",
                repeat="GlyphCount",
                aux=0,
                description="Array of GlyphIDs-in numerical order",
            ),
        ],
    ),
    (
        "CoverageFormat2",
        [
            FieldSpec(
                "uint16", "CoverageFormat", description="Format identifier-format = 2"
            ),
            FieldSpec("uint16", "RangeCount", description="Number of RangeRecords"),
            FieldSpec(
                "struct",
                "RangeRecord",
                repeat="RangeCount",
                aux=0,
                description="Array of glyph ranges-ordered by Start GlyphID",
            ),
        ],
    ),
    (
        "RangeRecord",
        [
            FieldSpec("GlyphID", "Start", description="First GlyphID in the range"),
            FieldSpec("GlyphID", "End", description="Last GlyphID in the range"),
            FieldSpec(
                "uint16",
                "StartCoverageIndex",
                description="Coverage Index of first GlyphID in range",
            ),
        ],
    ),
    (
        "ClassDefFormat1",
        [
            FieldSpec(
                "uint16", "ClassFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "GlyphID",
                "StartGlyph",
                description="First GlyphID of the ClassValueArray",
            ),
            FieldSpec(
                "uint16", "GlyphCount", description="Size of the ClassValueArray"
            ),
            FieldSpec(
                "uint16",
                "ClassValueArray",
                repeat="GlyphCount",
                aux=0,
                description="Array of Class Values-one per GlyphID",
            ),
        ],
    ),
    (
        "ClassDefFormat2",
        [
            FieldSpec(
                "uint16", "ClassFormat", description="Format identifier-format = 2"
            ),
            FieldSpec(
                "uint16", "ClassRangeCount", description="Number of ClassRangeRecords"
            ),
            FieldSpec(
                "struct",
                "ClassRangeRecord",
                repeat="ClassRangeCount",
                aux=0,
                description="Array of ClassRangeRecords-ordered by Start GlyphID",
            ),
        ],
    ),
    (
        "ClassRangeRecord",
        [
            FieldSpec("GlyphID", "Start", description="First GlyphID in the range"),
            FieldSpec("GlyphID", "End", description="Last GlyphID in the range"),
            FieldSpec(
                "uint16", "Class", description="Applied to all glyphs in the range"
            ),
        ],
    ),
    (
        "Device",
        [
            FieldSpec(
                "uint16", "StartSize", description="Smallest size to correct-in ppem"
            ),
            FieldSpec(
                "uint16", "EndSize", description="Largest size to correct-in ppem"
            ),
            FieldSpec(
                "uint16",
                "DeltaFormat",
                description="Format of DeltaValue array data: 1, 2, or 3",
            ),
            FieldSpec(
                "DeltaValue",
                "DeltaValue",
                aux="DeltaFormat in (1,2,3)",
                description="Array of compressed data",
            ),
        ],
    ),
    #
    # gpos
    #
    (
        "GPOS",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version of the GPOS table- 0x00010000 or 0x00010001",
            ),
            FieldSpec(
                "Offset",
                "ScriptList",
                description="Offset to ScriptList table-from beginning of GPOS table",
            ),
            FieldSpec(
                "Offset",
                "FeatureList",
                description="Offset to FeatureList table-from beginning of GPOS table",
            ),
            FieldSpec(
                "Offset",
                "LookupList",
                description="Offset to LookupList table-from beginning of GPOS table",
            ),
            FieldSpec(
                "LOffset",
                "FeatureVariations",
                aux="Version >= 0x00010001",
                description="Offset to FeatureVariations table-from beginning of GPOS table",
            ),
        ],
    ),
    (
        "SinglePosFormat1",
        [
            FieldSpec(
                "uint16", "PosFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of SinglePos subtable",
            ),
            FieldSpec(
                "uint16",
                "ValueFormat",
                description="Defines the types of data in the ValueRecord",
            ),
            FieldSpec(
                "ValueRecord",
                "Value",
                description="Defines positioning value(s)-applied to all glyphs in the Coverage table",
            ),
        ],
    ),
    (
        "SinglePosFormat2",
        [
            FieldSpec(
                "uint16", "PosFormat", description="Format identifier-format = 2"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of SinglePos subtable",
            ),
            FieldSpec(
                "uint16",
                "ValueFormat",
                description="Defines the types of data in the ValueRecord",
            ),
            FieldSpec("uint16", "ValueCount", description="Number of ValueRecords"),
            FieldSpec(
                "ValueRecord",
                "Value",
                repeat="ValueCount",
                aux=0,
                description="Array of ValueRecords-positioning values applied to glyphs",
            ),
        ],
    ),
    (
        "PairPosFormat1",
        [
            FieldSpec(
                "uint16", "PosFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of PairPos subtable-only the first glyph in each pair",
            ),
            FieldSpec(
                "uint16",
                "ValueFormat1",
                description="Defines the types of data in ValueRecord1-for the first glyph in the pair -may be zero (0)",
            ),
            FieldSpec(
                "uint16",
                "ValueFormat2",
                description="Defines the types of data in ValueRecord2-for the second glyph in the pair -may be zero (0)",
            ),
            FieldSpec("uint16", "PairSetCount", description="Number of PairSet tables"),
            FieldSpec(
                "Offset",
                "PairSet",
                repeat="PairSetCount",
                aux=0,
                description="Array of offsets to PairSet tables-from beginning of PairPos subtable-ordered by Coverage Index",
            ),
        ],
    ),
    (
        "PairSet",
        [
            FieldSpec(
                "uint16", "PairValueCount", description="Number of PairValueRecords"
            ),
            FieldSpec(
                "struct",
                "PairValueRecord",
                repeat="PairValueCount",
                aux=0,
                description="Array of PairValueRecords-ordered by GlyphID of the second glyph",
            ),
        ],
    ),
    (
        "PairValueRecord",
        [
            FieldSpec(
                "GlyphID",
                "SecondGlyph",
                description="GlyphID of second glyph in the pair-first glyph is listed in the Coverage table",
            ),
            FieldSpec(
                "ValueRecord",
                "Value1",
                description="Positioning data for the first glyph in the pair",
            ),
            FieldSpec(
                "ValueRecord",
                "Value2",
                description="Positioning data for the second glyph in the pair",
            ),
        ],
    ),
    (
        "PairPosFormat2",
        [
            FieldSpec(
                "uint16", "PosFormat", description="Format identifier-format = 2"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of PairPos subtable-for the first glyph of the pair",
            ),
            FieldSpec(
                "uint16",
                "ValueFormat1",
                description="ValueRecord definition-for the first glyph of the pair-may be zero (0)",
            ),
            FieldSpec(
                "uint16",
                "ValueFormat2",
                description="ValueRecord definition-for the second glyph of the pair-may be zero (0)",
            ),
            FieldSpec(
                "Offset",
                "ClassDef1",
                description="Offset to ClassDef table-from beginning of PairPos subtable-for the first glyph of the pair",
            ),
            FieldSpec(
                "Offset",
                "ClassDef2",
                description="Offset to ClassDef table-from beginning of PairPos subtable-for the second glyph of the pair",
            ),
            FieldSpec(
                "uint16",
                "Class1Count",
                description="Number of classes in ClassDef1 table-includes Class0",
            ),
            FieldSpec(
                "uint16",
                "Class2Count",
                description="Number of classes in ClassDef2 table-includes Class0",
            ),
            FieldSpec(
                "struct",
                "Class1Record",
                repeat="Class1Count",
                aux=0,
                description="Array of Class1 records-ordered by Class1",
            ),
        ],
    ),
    (
        "Class1Record",
        [
            FieldSpec(
                "struct",
                "Class2Record",
                repeat="Class2Count",
                aux=0,
                description="Array of Class2 records-ordered by Class2",
            ),
        ],
    ),
    (
        "Class2Record",
        [
            FieldSpec(
                "ValueRecord",
                "Value1",
                description="Positioning for first glyph-empty if ValueFormat1 = 0",
            ),
            FieldSpec(
                "ValueRecord",
                "Value2",
                description="Positioning for second glyph-empty if ValueFormat2 = 0",
            ),
        ],
    ),
    (
        "CursivePosFormat1",
        [
            FieldSpec(
                "uint16", "PosFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of CursivePos subtable",
            ),
            FieldSpec(
                "uint16", "EntryExitCount", description="Number of EntryExit records"
            ),
            FieldSpec(
                "struct",
                "EntryExitRecord",
                repeat="EntryExitCount",
                aux=0,
                description="Array of EntryExit records-in Coverage Index order",
            ),
        ],
    ),
    (
        "EntryExitRecord",
        [
            FieldSpec(
                "Offset",
                "EntryAnchor",
                description="Offset to EntryAnchor table-from beginning of CursivePos subtable-may be NULL",
            ),
            FieldSpec(
                "Offset",
                "ExitAnchor",
                description="Offset to ExitAnchor table-from beginning of CursivePos subtable-may be NULL",
            ),
        ],
    ),
    (
        "MarkBasePosFormat1",
        [
            FieldSpec(
                "uint16", "PosFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "Offset",
                "MarkCoverage",
                description="Offset to MarkCoverage table-from beginning of MarkBasePos subtable",
            ),
            FieldSpec(
                "Offset",
                "BaseCoverage",
                description="Offset to BaseCoverage table-from beginning of MarkBasePos subtable",
            ),
            FieldSpec(
                "uint16",
                "ClassCount",
                description="Number of classes defined for marks",
            ),
            FieldSpec(
                "Offset",
                "MarkArray",
                description="Offset to MarkArray table-from beginning of MarkBasePos subtable",
            ),
            FieldSpec(
                "Offset",
                "BaseArray",
                description="Offset to BaseArray table-from beginning of MarkBasePos subtable",
            ),
        ],
    ),
    (
        "BaseArray",
        [
            FieldSpec("uint16", "BaseCount", description="Number of BaseRecords"),
            FieldSpec(
                "struct",
                "BaseRecord",
                repeat="BaseCount",
                aux=0,
                description="Array of BaseRecords-in order of BaseCoverage Index",
            ),
        ],
    ),
    (
        "BaseRecord",
        [
            FieldSpec(
                "Offset",
                "BaseAnchor",
                repeat="ClassCount",
                aux=0,
                description="Array of offsets (one per class) to Anchor tables-from beginning of BaseArray table-ordered by class-zero-based",
            ),
        ],
    ),
    (
        "MarkLigPosFormat1",
        [
            FieldSpec(
                "uint16", "PosFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "Offset",
                "MarkCoverage",
                description="Offset to Mark Coverage table-from beginning of MarkLigPos subtable",
            ),
            FieldSpec(
                "Offset",
                "LigatureCoverage",
                description="Offset to Ligature Coverage table-from beginning of MarkLigPos subtable",
            ),
            FieldSpec(
                "uint16", "ClassCount", description="Number of defined mark classes"
            ),
            FieldSpec(
                "Offset",
                "MarkArray",
                description="Offset to MarkArray table-from beginning of MarkLigPos subtable",
            ),
            FieldSpec(
                "Offset",
                "LigatureArray",
                description="Offset to LigatureArray table-from beginning of MarkLigPos subtable",
            ),
        ],
    ),
    (
        "LigatureArray",
        [
            FieldSpec(
                "uint16",
                "LigatureCount",
                description="Number of LigatureAttach table offsets",
            ),
            FieldSpec(
                "Offset",
                "LigatureAttach",
                repeat="LigatureCount",
                aux=0,
                description="Array of offsets to LigatureAttach tables-from beginning of LigatureArray table-ordered by LigatureCoverage Index",
            ),
        ],
    ),
    (
        "LigatureAttach",
        [
            FieldSpec(
                "uint16",
                "ComponentCount",
                description="Number of ComponentRecords in this ligature",
            ),
            FieldSpec(
                "struct",
                "ComponentRecord",
                repeat="ComponentCount",
                aux=0,
                description="Array of Component records-ordered in writing direction",
            ),
        ],
    ),
    (
        "ComponentRecord",
        [
            FieldSpec(
                "Offset",
                "LigatureAnchor",
                repeat="ClassCount",
                aux=0,
                description="Array of offsets (one per class) to Anchor tables-from beginning of LigatureAttach table-ordered by class-NULL if a component does not have an attachment for a class-zero-based array",
            ),
        ],
    ),
    (
        "MarkMarkPosFormat1",
        [
            FieldSpec(
                "uint16", "PosFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "Offset",
                "Mark1Coverage",
                description="Offset to Combining Mark Coverage table-from beginning of MarkMarkPos subtable",
            ),
            FieldSpec(
                "Offset",
                "Mark2Coverage",
                description="Offset to Base Mark Coverage table-from beginning of MarkMarkPos subtable",
            ),
            FieldSpec(
                "uint16",
                "ClassCount",
                description="Number of Combining Mark classes defined",
            ),
            FieldSpec(
                "Offset",
                "Mark1Array",
                description="Offset to MarkArray table for Mark1-from beginning of MarkMarkPos subtable",
            ),
            FieldSpec(
                "Offset",
                "Mark2Array",
                description="Offset to Mark2Array table for Mark2-from beginning of MarkMarkPos subtable",
            ),
        ],
    ),
    (
        "Mark2Array",
        [
            FieldSpec("uint16", "Mark2Count", description="Number of Mark2 records"),
            FieldSpec(
                "struct",
                "Mark2Record",
                repeat="Mark2Count",
                aux=0,
                description="Array of Mark2 records-in Coverage order",
            ),
        ],
    ),
    (
        "Mark2Record",
        [
            FieldSpec(
                "Offset",
                "Mark2Anchor",
                repeat="ClassCount",
                aux=0,
                description="Array of offsets (one per class) to Anchor tables-from beginning of Mark2Array table-zero-based array",
            ),
        ],
    ),
    (
        "PosLookupRecord",
        [
            FieldSpec(
                "uint16",
                "SequenceIndex",
                description="Index to input glyph sequence-first glyph = 0",
            ),
            FieldSpec(
                "uint16",
                "LookupListIndex",
                description="Lookup to apply to that position-zero-based",
            ),
        ],
    ),
    (
        "ContextPosFormat1",
        [
            FieldSpec(
                "uint16", "PosFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of ContextPos subtable",
            ),
            FieldSpec(
                "uint16", "PosRuleSetCount", description="Number of PosRuleSet tables"
            ),
            FieldSpec(
                "Offset",
                "PosRuleSet",
                repeat="PosRuleSetCount",
                aux=0,
                description="Array of offsets to PosRuleSet tables-from beginning of ContextPos subtable-ordered by Coverage Index",
            ),
        ],
    ),
    (
        "PosRuleSet",
        [
            FieldSpec("uint16", "PosRuleCount", description="Number of PosRule tables"),
            FieldSpec(
                "Offset",
                "PosRule",
                repeat="PosRuleCount",
                aux=0,
                description="Array of offsets to PosRule tables-from beginning of PosRuleSet-ordered by preference",
            ),
        ],
    ),
    (
        "PosRule",
        [
            FieldSpec(
                "uint16",
                "GlyphCount",
                description="Number of glyphs in the Input glyph sequence",
            ),
            FieldSpec("uint16", "PosCount", description="Number of PosLookupRecords"),
            FieldSpec(
                "GlyphID",
                "Input",
                repeat="GlyphCount",
                aux=-1,
                description="Array of input GlyphIDs-starting with the second glyph",
            ),
            FieldSpec(
                "struct",
                "PosLookupRecord",
                repeat="PosCount",
                aux=0,
                description="Array of positioning lookups-in design order",
            ),
        ],
    ),
    (
        "ContextPosFormat2",
        [
            FieldSpec(
                "uint16", "PosFormat", description="Format identifier-format = 2"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of ContextPos subtable",
            ),
            FieldSpec(
                "Offset",
                "ClassDef",
                description="Offset to ClassDef table-from beginning of ContextPos subtable",
            ),
            FieldSpec(
                "uint16", "PosClassSetCount", description="Number of PosClassSet tables"
            ),
            FieldSpec(
                "Offset",
                "PosClassSet",
                repeat="PosClassSetCount",
                aux=0,
                description="Array of offsets to PosClassSet tables-from beginning of ContextPos subtable-ordered by class-may be NULL",
            ),
        ],
    ),
    (
        "PosClassSet",
        [
            FieldSpec(
                "uint16",
                "PosClassRuleCount",
                description="Number of PosClassRule tables",
            ),
            FieldSpec(
                "Offset",
                "PosClassRule",
                repeat="PosClassRuleCount",
                aux=0,
                description="Array of offsets to PosClassRule tables-from beginning of PosClassSet-ordered by preference",
            ),
        ],
    ),
    (
        "PosClassRule",
        [
            FieldSpec(
                "uint16", "GlyphCount", description="Number of glyphs to be matched"
            ),
            FieldSpec("uint16", "PosCount", description="Number of PosLookupRecords"),
            FieldSpec(
                "uint16",
                "Class",
                repeat="GlyphCount",
                aux=-1,
                description="Array of classes-beginning with the second class-to be matched to the input glyph sequence",
            ),
            FieldSpec(
                "struct",
                "PosLookupRecord",
                repeat="PosCount",
                aux=0,
                description="Array of positioning lookups-in design order",
            ),
        ],
    ),
    (
        "ContextPosFormat3",
        [
            FieldSpec(
                "uint16", "PosFormat", description="Format identifier-format = 3"
            ),
            FieldSpec(
                "uint16",
                "GlyphCount",
                description="Number of glyphs in the input sequence",
            ),
            FieldSpec("uint16", "PosCount", description="Number of PosLookupRecords"),
            FieldSpec(
                "Offset",
                "Coverage",
                repeat="GlyphCount",
                aux=0,
                description="Array of offsets to Coverage tables-from beginning of ContextPos subtable",
            ),
            FieldSpec(
                "struct",
                "PosLookupRecord",
                repeat="PosCount",
                aux=0,
                description="Array of positioning lookups-in design order",
            ),
        ],
    ),
    (
        "ChainContextPosFormat1",
        [
            FieldSpec(
                "uint16", "PosFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of ContextPos subtable",
            ),
            FieldSpec(
                "uint16",
                "ChainPosRuleSetCount",
                description="Number of ChainPosRuleSet tables",
            ),
            FieldSpec(
                "Offset",
                "ChainPosRuleSet",
                repeat="ChainPosRuleSetCount",
                aux=0,
                description="Array of offsets to ChainPosRuleSet tables-from beginning of ContextPos subtable-ordered by Coverage Index",
            ),
        ],
    ),
    (
        "ChainPosRuleSet",
        [
            FieldSpec(
                "uint16",
                "ChainPosRuleCount",
                description="Number of ChainPosRule tables",
            ),
            FieldSpec(
                "Offset",
                "ChainPosRule",
                repeat="ChainPosRuleCount",
                aux=0,
                description="Array of offsets to ChainPosRule tables-from beginning of ChainPosRuleSet-ordered by preference",
            ),
        ],
    ),
    (
        "ChainPosRule",
        [
            FieldSpec(
                "uint16",
                "BacktrackGlyphCount",
                description="Total number of glyphs in the backtrack sequence (number of glyphs to be matched before the first glyph)",
            ),
            FieldSpec(
                "GlyphID",
                "Backtrack",
                repeat="BacktrackGlyphCount",
                aux=0,
                description="Array of backtracking GlyphID's (to be matched before the input sequence)",
            ),
            FieldSpec(
                "uint16",
                "InputGlyphCount",
                description="Total number of glyphs in the input sequence (includes the first glyph)",
            ),
            FieldSpec(
                "GlyphID",
                "Input",
                repeat="InputGlyphCount",
                aux=-1,
                description="Array of input GlyphIDs (start with second glyph)",
            ),
            FieldSpec(
                "uint16",
                "LookAheadGlyphCount",
                description="Total number of glyphs in the look ahead sequence (number of glyphs to be matched after the input sequence)",
            ),
            FieldSpec(
                "GlyphID",
                "LookAhead",
                repeat="LookAheadGlyphCount",
                aux=0,
                description="Array of lookahead GlyphID's (to be matched after the input sequence)",
            ),
            FieldSpec("uint16", "PosCount", description="Number of PosLookupRecords"),
            FieldSpec(
                "struct",
                "PosLookupRecord",
                repeat="PosCount",
                aux=0,
                description="Array of PosLookupRecords (in design order)",
            ),
        ],
    ),
    (
        "ChainContextPosFormat2",
        [
            FieldSpec(
                "uint16", "PosFormat", description="Format identifier-format = 2"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of ChainContextPos subtable",
            ),
            FieldSpec(
                "Offset",
                "BacktrackClassDef",
                description="Offset to ClassDef table containing backtrack sequence context-from beginning of ChainContextPos subtable",
            ),
            FieldSpec(
                "Offset",
                "InputClassDef",
                description="Offset to ClassDef table containing input sequence context-from beginning of ChainContextPos subtable",
            ),
            FieldSpec(
                "Offset",
                "LookAheadClassDef",
                description="Offset to ClassDef table containing lookahead sequence context-from beginning of ChainContextPos subtable",
            ),
            FieldSpec(
                "uint16",
                "ChainPosClassSetCount",
                description="Number of ChainPosClassSet tables",
            ),
            FieldSpec(
                "Offset",
                "ChainPosClassSet",
                repeat="ChainPosClassSetCount",
                aux=0,
                description="Array of offsets to ChainPosClassSet tables-from beginning of ChainContextPos subtable-ordered by input class-may be NULL",
            ),
        ],
    ),
    (
        "ChainPosClassSet",
        [
            FieldSpec(
                "uint16",
                "ChainPosClassRuleCount",
                description="Number of ChainPosClassRule tables",
            ),
            FieldSpec(
                "Offset",
                "ChainPosClassRule",
                repeat="ChainPosClassRuleCount",
                aux=0,
                description="Array of offsets to ChainPosClassRule tables-from beginning of ChainPosClassSet-ordered by preference",
            ),
        ],
    ),
    (
        "ChainPosClassRule",
        [
            FieldSpec(
                "uint16",
                "BacktrackGlyphCount",
                description="Total number of glyphs in the backtrack sequence (number of glyphs to be matched before the first glyph)",
            ),
            FieldSpec(
                "uint16",
                "Backtrack",
                repeat="BacktrackGlyphCount",
                aux=0,
                description="Array of backtracking classes(to be matched before the input sequence)",
            ),
            FieldSpec(
                "uint16",
                "InputGlyphCount",
                description="Total number of classes in the input sequence (includes the first class)",
            ),
            FieldSpec(
                "uint16",
                "Input",
                repeat="InputGlyphCount",
                aux=-1,
                description="Array of input classes(start with second class; to be matched with the input glyph sequence)",
            ),
            FieldSpec(
                "uint16",
                "LookAheadGlyphCount",
                description="Total number of classes in the look ahead sequence (number of classes to be matched after the input sequence)",
            ),
            FieldSpec(
                "uint16",
                "LookAhead",
                repeat="LookAheadGlyphCount",
                aux=0,
                description="Array of lookahead classes(to be matched after the input sequence)",
            ),
            FieldSpec("uint16", "PosCount", description="Number of PosLookupRecords"),
            FieldSpec(
                "struct",
                "PosLookupRecord",
                repeat="PosCount",
                aux=0,
                description="Array of PosLookupRecords (in design order)",
            ),
        ],
    ),
    (
        "ChainContextPosFormat3",
        [
            FieldSpec(
                "uint16", "PosFormat", description="Format identifier-format = 3"
            ),
            FieldSpec(
                "uint16",
                "BacktrackGlyphCount",
                description="Number of glyphs in the backtracking sequence",
            ),
            FieldSpec(
                "Offset",
                "BacktrackCoverage",
                repeat="BacktrackGlyphCount",
                aux=0,
                description="Array of offsets to coverage tables in backtracking sequence, in glyph sequence order",
            ),
            FieldSpec(
                "uint16",
                "InputGlyphCount",
                description="Number of glyphs in input sequence",
            ),
            FieldSpec(
                "Offset",
                "InputCoverage",
                repeat="InputGlyphCount",
                aux=0,
                description="Array of offsets to coverage tables in input sequence, in glyph sequence order",
            ),
            FieldSpec(
                "uint16",
                "LookAheadGlyphCount",
                description="Number of glyphs in lookahead sequence",
            ),
            FieldSpec(
                "Offset",
                "LookAheadCoverage",
                repeat="LookAheadGlyphCount",
                aux=0,
                description="Array of offsets to coverage tables in lookahead sequence, in glyph sequence order",
            ),
            FieldSpec("uint16", "PosCount", description="Number of PosLookupRecords"),
            FieldSpec(
                "struct",
                "PosLookupRecord",
                repeat="PosCount",
                aux=0,
                description="Array of PosLookupRecords,in design order",
            ),
        ],
    ),
    (
        "ExtensionPosFormat1",
        [
            FieldSpec(
                "uint16", "ExtFormat", description="Format identifier. Set to 1."
            ),
            FieldSpec(
                "uint16",
                "ExtensionLookupType",
                description="Lookup type of subtable referenced by ExtensionOffset (i.e. the extension subtable).",
            ),
            FieldSpec("LOffset", "ExtSubTable", description="Offset to SubTable"),
        ],
    ),
    # 	('ValueRecord', [
    # 		('int16', 'XPlacement', None, None, 'Horizontal adjustment for placement-in design units'),
    # 		('int16', 'YPlacement', None, None, 'Vertical adjustment for placement-in design units'),
    # 		('int16', 'XAdvance', None, None, 'Horizontal adjustment for advance-in design units (only used for horizontal writing)'),
    # 		('int16', 'YAdvance', None, None, 'Vertical adjustment for advance-in design units (only used for vertical writing)'),
    # 		('Offset', 'XPlaDevice', None, None, 'Offset to Device table for horizontal placement-measured from beginning of PosTable (may be NULL)'),
    # 		('Offset', 'YPlaDevice', None, None, 'Offset to Device table for vertical placement-measured from beginning of PosTable (may be NULL)'),
    # 		('Offset', 'XAdvDevice', None, None, 'Offset to Device table for horizontal advance-measured from beginning of PosTable (may be NULL)'),
    # 		('Offset', 'YAdvDevice', None, None, 'Offset to Device table for vertical advance-measured from beginning of PosTable (may be NULL)'),
    # 	]),
    (
        "AnchorFormat1",
        [
            FieldSpec(
                "uint16", "AnchorFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "int16", "XCoordinate", description="Horizontal value-in design units"
            ),
            FieldSpec(
                "int16", "YCoordinate", description="Vertical value-in design units"
            ),
        ],
    ),
    (
        "AnchorFormat2",
        [
            FieldSpec(
                "uint16", "AnchorFormat", description="Format identifier-format = 2"
            ),
            FieldSpec(
                "int16", "XCoordinate", description="Horizontal value-in design units"
            ),
            FieldSpec(
                "int16", "YCoordinate", description="Vertical value-in design units"
            ),
            FieldSpec(
                "uint16", "AnchorPoint", description="Index to glyph contour point"
            ),
        ],
    ),
    (
        "AnchorFormat3",
        [
            FieldSpec(
                "uint16", "AnchorFormat", description="Format identifier-format = 3"
            ),
            FieldSpec(
                "int16", "XCoordinate", description="Horizontal value-in design units"
            ),
            FieldSpec(
                "int16", "YCoordinate", description="Vertical value-in design units"
            ),
            FieldSpec(
                "Offset",
                "XDeviceTable",
                description="Offset to Device table for X coordinate- from beginning of Anchor table (may be NULL)",
            ),
            FieldSpec(
                "Offset",
                "YDeviceTable",
                description="Offset to Device table for Y coordinate- from beginning of Anchor table (may be NULL)",
            ),
        ],
    ),
    (
        "MarkArray",
        [
            FieldSpec("uint16", "MarkCount", description="Number of MarkRecords"),
            FieldSpec(
                "struct",
                "MarkRecord",
                repeat="MarkCount",
                aux=0,
                description="Array of MarkRecords-in Coverage order",
            ),
        ],
    ),
    (
        "MarkRecord",
        [
            FieldSpec("uint16", "Class", description="Class defined for this mark"),
            FieldSpec(
                "Offset",
                "MarkAnchor",
                description="Offset to Anchor table-from beginning of MarkArray table",
            ),
        ],
    ),
    #
    # gsub
    #
    (
        "GSUB",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version of the GSUB table- 0x00010000 or 0x00010001",
            ),
            FieldSpec(
                "Offset",
                "ScriptList",
                description="Offset to ScriptList table-from beginning of GSUB table",
            ),
            FieldSpec(
                "Offset",
                "FeatureList",
                description="Offset to FeatureList table-from beginning of GSUB table",
            ),
            FieldSpec(
                "Offset",
                "LookupList",
                description="Offset to LookupList table-from beginning of GSUB table",
            ),
            FieldSpec(
                "LOffset",
                "FeatureVariations",
                aux="Version >= 0x00010001",
                description="Offset to FeatureVariations table-from beginning of GSUB table",
            ),
        ],
    ),
    (
        "SingleSubstFormat1",
        [
            FieldSpec(
                "uint16", "SubstFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of Substitution table",
            ),
            FieldSpec(
                "uint16",
                "DeltaGlyphID",
                description="Add to original GlyphID modulo 65536 to get substitute GlyphID",
            ),
        ],
    ),
    (
        "SingleSubstFormat2",
        [
            FieldSpec(
                "uint16", "SubstFormat", description="Format identifier-format = 2"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of Substitution table",
            ),
            FieldSpec(
                "uint16",
                "GlyphCount",
                description="Number of GlyphIDs in the Substitute array",
            ),
            FieldSpec(
                "GlyphID",
                "Substitute",
                repeat="GlyphCount",
                aux=0,
                description="Array of substitute GlyphIDs-ordered by Coverage Index",
            ),
        ],
    ),
    (
        "MultipleSubstFormat1",
        [
            FieldSpec(
                "uint16", "SubstFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of Substitution table",
            ),
            FieldSpec(
                "uint16",
                "SequenceCount",
                description="Number of Sequence table offsets in the Sequence array",
            ),
            FieldSpec(
                "Offset",
                "Sequence",
                repeat="SequenceCount",
                aux=0,
                description="Array of offsets to Sequence tables-from beginning of Substitution table-ordered by Coverage Index",
            ),
        ],
    ),
    (
        "Sequence",
        [
            FieldSpec(
                "uint16",
                "GlyphCount",
                description="Number of GlyphIDs in the Substitute array. This should always be greater than 0.",
            ),
            FieldSpec(
                "GlyphID",
                "Substitute",
                repeat="GlyphCount",
                aux=0,
                description="String of GlyphIDs to substitute",
            ),
        ],
    ),
    (
        "AlternateSubstFormat1",
        [
            FieldSpec(
                "uint16", "SubstFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of Substitution table",
            ),
            FieldSpec(
                "uint16",
                "AlternateSetCount",
                description="Number of AlternateSet tables",
            ),
            FieldSpec(
                "Offset",
                "AlternateSet",
                repeat="AlternateSetCount",
                aux=0,
                description="Array of offsets to AlternateSet tables-from beginning of Substitution table-ordered by Coverage Index",
            ),
        ],
    ),
    (
        "AlternateSet",
        [
            FieldSpec(
                "uint16",
                "GlyphCount",
                description="Number of GlyphIDs in the Alternate array",
            ),
            FieldSpec(
                "GlyphID",
                "Alternate",
                repeat="GlyphCount",
                aux=0,
                description="Array of alternate GlyphIDs-in arbitrary order",
            ),
        ],
    ),
    (
        "LigatureSubstFormat1",
        [
            FieldSpec(
                "uint16", "SubstFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of Substitution table",
            ),
            FieldSpec(
                "uint16", "LigSetCount", description="Number of LigatureSet tables"
            ),
            FieldSpec(
                "Offset",
                "LigatureSet",
                repeat="LigSetCount",
                aux=0,
                description="Array of offsets to LigatureSet tables-from beginning of Substitution table-ordered by Coverage Index",
            ),
        ],
    ),
    (
        "LigatureSet",
        [
            FieldSpec(
                "uint16", "LigatureCount", description="Number of Ligature tables"
            ),
            FieldSpec(
                "Offset",
                "Ligature",
                repeat="LigatureCount",
                aux=0,
                description="Array of offsets to Ligature tables-from beginning of LigatureSet table-ordered by preference",
            ),
        ],
    ),
    (
        "Ligature",
        [
            FieldSpec(
                "GlyphID", "LigGlyph", description="GlyphID of ligature to substitute"
            ),
            FieldSpec(
                "uint16",
                "CompCount",
                description="Number of components in the ligature",
            ),
            FieldSpec(
                "GlyphID",
                "Component",
                repeat="CompCount",
                aux=-1,
                description="Array of component GlyphIDs-start with the second component-ordered in writing direction",
            ),
        ],
    ),
    (
        "SubstLookupRecord",
        [
            FieldSpec(
                "uint16",
                "SequenceIndex",
                description="Index into current glyph sequence-first glyph = 0",
            ),
            FieldSpec(
                "uint16",
                "LookupListIndex",
                description="Lookup to apply to that position-zero-based",
            ),
        ],
    ),
    (
        "ContextSubstFormat1",
        [
            FieldSpec(
                "uint16", "SubstFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of Substitution table",
            ),
            FieldSpec(
                "uint16",
                "SubRuleSetCount",
                description="Number of SubRuleSet tables-must equal GlyphCount in Coverage table",
            ),
            FieldSpec(
                "Offset",
                "SubRuleSet",
                repeat="SubRuleSetCount",
                aux=0,
                description="Array of offsets to SubRuleSet tables-from beginning of Substitution table-ordered by Coverage Index",
            ),
        ],
    ),
    (
        "SubRuleSet",
        [
            FieldSpec("uint16", "SubRuleCount", description="Number of SubRule tables"),
            FieldSpec(
                "Offset",
                "SubRule",
                repeat="SubRuleCount",
                aux=0,
                description="Array of offsets to SubRule tables-from beginning of SubRuleSet table-ordered by preference",
            ),
        ],
    ),
    (
        "SubRule",
        [
            FieldSpec(
                "uint16",
                "GlyphCount",
                description="Total number of glyphs in input glyph sequence-includes the first glyph",
            ),
            FieldSpec(
                "uint16", "SubstCount", description="Number of SubstLookupRecords"
            ),
            FieldSpec(
                "GlyphID",
                "Input",
                repeat="GlyphCount",
                aux=-1,
                description="Array of input GlyphIDs-start with second glyph",
            ),
            FieldSpec(
                "struct",
                "SubstLookupRecord",
                repeat="SubstCount",
                aux=0,
                description="Array of SubstLookupRecords-in design order",
            ),
        ],
    ),
    (
        "ContextSubstFormat2",
        [
            FieldSpec(
                "uint16", "SubstFormat", description="Format identifier-format = 2"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of Substitution table",
            ),
            FieldSpec(
                "Offset",
                "ClassDef",
                description="Offset to glyph ClassDef table-from beginning of Substitution table",
            ),
            FieldSpec(
                "uint16", "SubClassSetCount", description="Number of SubClassSet tables"
            ),
            FieldSpec(
                "Offset",
                "SubClassSet",
                repeat="SubClassSetCount",
                aux=0,
                description="Array of offsets to SubClassSet tables-from beginning of Substitution table-ordered by class-may be NULL",
            ),
        ],
    ),
    (
        "SubClassSet",
        [
            FieldSpec(
                "uint16",
                "SubClassRuleCount",
                description="Number of SubClassRule tables",
            ),
            FieldSpec(
                "Offset",
                "SubClassRule",
                repeat="SubClassRuleCount",
                aux=0,
                description="Array of offsets to SubClassRule tables-from beginning of SubClassSet-ordered by preference",
            ),
        ],
    ),
    (
        "SubClassRule",
        [
            FieldSpec(
                "uint16",
                "GlyphCount",
                description="Total number of classes specified for the context in the rule-includes the first class",
            ),
            FieldSpec(
                "uint16", "SubstCount", description="Number of SubstLookupRecords"
            ),
            FieldSpec(
                "uint16",
                "Class",
                repeat="GlyphCount",
                aux=-1,
                description="Array of classes-beginning with the second class-to be matched to the input glyph class sequence",
            ),
            FieldSpec(
                "struct",
                "SubstLookupRecord",
                repeat="SubstCount",
                aux=0,
                description="Array of Substitution lookups-in design order",
            ),
        ],
    ),
    (
        "ContextSubstFormat3",
        [
            FieldSpec(
                "uint16", "SubstFormat", description="Format identifier-format = 3"
            ),
            FieldSpec(
                "uint16",
                "GlyphCount",
                description="Number of glyphs in the input glyph sequence",
            ),
            FieldSpec(
                "uint16", "SubstCount", description="Number of SubstLookupRecords"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                repeat="GlyphCount",
                aux=0,
                description="Array of offsets to Coverage table-from beginning of Substitution table-in glyph sequence order",
            ),
            FieldSpec(
                "struct",
                "SubstLookupRecord",
                repeat="SubstCount",
                aux=0,
                description="Array of SubstLookupRecords-in design order",
            ),
        ],
    ),
    (
        "ChainContextSubstFormat1",
        [
            FieldSpec(
                "uint16", "SubstFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of Substitution table",
            ),
            FieldSpec(
                "uint16",
                "ChainSubRuleSetCount",
                description="Number of ChainSubRuleSet tables-must equal GlyphCount in Coverage table",
            ),
            FieldSpec(
                "Offset",
                "ChainSubRuleSet",
                repeat="ChainSubRuleSetCount",
                aux=0,
                description="Array of offsets to ChainSubRuleSet tables-from beginning of Substitution table-ordered by Coverage Index",
            ),
        ],
    ),
    (
        "ChainSubRuleSet",
        [
            FieldSpec(
                "uint16",
                "ChainSubRuleCount",
                description="Number of ChainSubRule tables",
            ),
            FieldSpec(
                "Offset",
                "ChainSubRule",
                repeat="ChainSubRuleCount",
                aux=0,
                description="Array of offsets to ChainSubRule tables-from beginning of ChainSubRuleSet table-ordered by preference",
            ),
        ],
    ),
    (
        "ChainSubRule",
        [
            FieldSpec(
                "uint16",
                "BacktrackGlyphCount",
                description="Total number of glyphs in the backtrack sequence (number of glyphs to be matched before the first glyph)",
            ),
            FieldSpec(
                "GlyphID",
                "Backtrack",
                repeat="BacktrackGlyphCount",
                aux=0,
                description="Array of backtracking GlyphID's (to be matched before the input sequence)",
            ),
            FieldSpec(
                "uint16",
                "InputGlyphCount",
                description="Total number of glyphs in the input sequence (includes the first glyph)",
            ),
            FieldSpec(
                "GlyphID",
                "Input",
                repeat="InputGlyphCount",
                aux=-1,
                description="Array of input GlyphIDs (start with second glyph)",
            ),
            FieldSpec(
                "uint16",
                "LookAheadGlyphCount",
                description="Total number of glyphs in the look ahead sequence (number of glyphs to be matched after the input sequence)",
            ),
            FieldSpec(
                "GlyphID",
                "LookAhead",
                repeat="LookAheadGlyphCount",
                aux=0,
                description="Array of lookahead GlyphID's (to be matched after the input sequence)",
            ),
            FieldSpec(
                "uint16", "SubstCount", description="Number of SubstLookupRecords"
            ),
            FieldSpec(
                "struct",
                "SubstLookupRecord",
                repeat="SubstCount",
                aux=0,
                description="Array of SubstLookupRecords (in design order)",
            ),
        ],
    ),
    (
        "ChainContextSubstFormat2",
        [
            FieldSpec(
                "uint16", "SubstFormat", description="Format identifier-format = 2"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table-from beginning of Substitution table",
            ),
            FieldSpec(
                "Offset",
                "BacktrackClassDef",
                description="Offset to glyph ClassDef table containing backtrack sequence data-from beginning of Substitution table",
            ),
            FieldSpec(
                "Offset",
                "InputClassDef",
                description="Offset to glyph ClassDef table containing input sequence data-from beginning of Substitution table",
            ),
            FieldSpec(
                "Offset",
                "LookAheadClassDef",
                description="Offset to glyph ClassDef table containing lookahead sequence data-from beginning of Substitution table",
            ),
            FieldSpec(
                "uint16",
                "ChainSubClassSetCount",
                description="Number of ChainSubClassSet tables",
            ),
            FieldSpec(
                "Offset",
                "ChainSubClassSet",
                repeat="ChainSubClassSetCount",
                aux=0,
                description="Array of offsets to ChainSubClassSet tables-from beginning of Substitution table-ordered by input class-may be NULL",
            ),
        ],
    ),
    (
        "ChainSubClassSet",
        [
            FieldSpec(
                "uint16",
                "ChainSubClassRuleCount",
                description="Number of ChainSubClassRule tables",
            ),
            FieldSpec(
                "Offset",
                "ChainSubClassRule",
                repeat="ChainSubClassRuleCount",
                aux=0,
                description="Array of offsets to ChainSubClassRule tables-from beginning of ChainSubClassSet-ordered by preference",
            ),
        ],
    ),
    (
        "ChainSubClassRule",
        [
            FieldSpec(
                "uint16",
                "BacktrackGlyphCount",
                description="Total number of glyphs in the backtrack sequence (number of glyphs to be matched before the first glyph)",
            ),
            FieldSpec(
                "uint16",
                "Backtrack",
                repeat="BacktrackGlyphCount",
                aux=0,
                description="Array of backtracking classes(to be matched before the input sequence)",
            ),
            FieldSpec(
                "uint16",
                "InputGlyphCount",
                description="Total number of classes in the input sequence (includes the first class)",
            ),
            FieldSpec(
                "uint16",
                "Input",
                repeat="InputGlyphCount",
                aux=-1,
                description="Array of input classes(start with second class; to be matched with the input glyph sequence)",
            ),
            FieldSpec(
                "uint16",
                "LookAheadGlyphCount",
                description="Total number of classes in the look ahead sequence (number of classes to be matched after the input sequence)",
            ),
            FieldSpec(
                "uint16",
                "LookAhead",
                repeat="LookAheadGlyphCount",
                aux=0,
                description="Array of lookahead classes(to be matched after the input sequence)",
            ),
            FieldSpec(
                "uint16", "SubstCount", description="Number of SubstLookupRecords"
            ),
            FieldSpec(
                "struct",
                "SubstLookupRecord",
                repeat="SubstCount",
                aux=0,
                description="Array of SubstLookupRecords (in design order)",
            ),
        ],
    ),
    (
        "ChainContextSubstFormat3",
        [
            FieldSpec(
                "uint16", "SubstFormat", description="Format identifier-format = 3"
            ),
            FieldSpec(
                "uint16",
                "BacktrackGlyphCount",
                description="Number of glyphs in the backtracking sequence",
            ),
            FieldSpec(
                "Offset",
                "BacktrackCoverage",
                repeat="BacktrackGlyphCount",
                aux=0,
                description="Array of offsets to coverage tables in backtracking sequence, in glyph sequence order",
            ),
            FieldSpec(
                "uint16",
                "InputGlyphCount",
                description="Number of glyphs in input sequence",
            ),
            FieldSpec(
                "Offset",
                "InputCoverage",
                repeat="InputGlyphCount",
                aux=0,
                description="Array of offsets to coverage tables in input sequence, in glyph sequence order",
            ),
            FieldSpec(
                "uint16",
                "LookAheadGlyphCount",
                description="Number of glyphs in lookahead sequence",
            ),
            FieldSpec(
                "Offset",
                "LookAheadCoverage",
                repeat="LookAheadGlyphCount",
                aux=0,
                description="Array of offsets to coverage tables in lookahead sequence, in glyph sequence order",
            ),
            FieldSpec(
                "uint16", "SubstCount", description="Number of SubstLookupRecords"
            ),
            FieldSpec(
                "struct",
                "SubstLookupRecord",
                repeat="SubstCount",
                aux=0,
                description="Array of SubstLookupRecords, in design order",
            ),
        ],
    ),
    (
        "ExtensionSubstFormat1",
        [
            FieldSpec(
                "uint16", "ExtFormat", description="Format identifier. Set to 1."
            ),
            FieldSpec(
                "uint16",
                "ExtensionLookupType",
                description="Lookup type of subtable referenced by ExtensionOffset (i.e. the extension subtable).",
            ),
            FieldSpec(
                "LOffset",
                "ExtSubTable",
                description="Array of offsets to Lookup tables-from beginning of LookupList -zero based (first lookup is Lookup index = 0)",
            ),
        ],
    ),
    (
        "ReverseChainSingleSubstFormat1",
        [
            FieldSpec(
                "uint16", "SubstFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "Offset",
                "Coverage",
                aux=0,
                description="Offset to Coverage table - from beginning of Substitution table",
            ),
            FieldSpec(
                "uint16",
                "BacktrackGlyphCount",
                description="Number of glyphs in the backtracking sequence",
            ),
            FieldSpec(
                "Offset",
                "BacktrackCoverage",
                repeat="BacktrackGlyphCount",
                aux=0,
                description="Array of offsets to coverage tables in backtracking sequence, in glyph sequence order",
            ),
            FieldSpec(
                "uint16",
                "LookAheadGlyphCount",
                description="Number of glyphs in lookahead sequence",
            ),
            FieldSpec(
                "Offset",
                "LookAheadCoverage",
                repeat="LookAheadGlyphCount",
                aux=0,
                description="Array of offsets to coverage tables in lookahead sequence, in glyph sequence order",
            ),
            FieldSpec(
                "uint16",
                "GlyphCount",
                description="Number of GlyphIDs in the Substitute array",
            ),
            FieldSpec(
                "GlyphID",
                "Substitute",
                repeat="GlyphCount",
                aux=0,
                description="Array of substitute GlyphIDs-ordered by Coverage index",
            ),
        ],
    ),
    #
    # gdef
    #
    (
        "GDEF",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version of the GDEF table- 0x00010000, 0x00010002, or 0x00010003",
            ),
            FieldSpec(
                "Offset",
                "GlyphClassDef",
                description="Offset to class definition table for glyph type-from beginning of GDEF header (may be NULL)",
            ),
            FieldSpec(
                "Offset",
                "AttachList",
                description="Offset to list of glyphs with attachment points-from beginning of GDEF header (may be NULL)",
            ),
            FieldSpec(
                "Offset",
                "LigCaretList",
                description="Offset to list of positioning points for ligature carets-from beginning of GDEF header (may be NULL)",
            ),
            FieldSpec(
                "Offset",
                "MarkAttachClassDef",
                description="Offset to class definition table for mark attachment type-from beginning of GDEF header (may be NULL)",
            ),
            FieldSpec(
                "Offset",
                "MarkGlyphSetsDef",
                aux="Version >= 0x00010002",
                description="Offset to the table of mark set definitions-from beginning of GDEF header (may be NULL)",
            ),
            FieldSpec(
                "LOffset",
                "VarStore",
                aux="Version >= 0x00010003",
                description="Offset to variation store (may be NULL)",
            ),
        ],
    ),
    (
        "AttachList",
        [
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table - from beginning of AttachList table",
            ),
            FieldSpec(
                "uint16",
                "GlyphCount",
                description="Number of glyphs with attachment points",
            ),
            FieldSpec(
                "Offset",
                "AttachPoint",
                repeat="GlyphCount",
                aux=0,
                description="Array of offsets to AttachPoint tables-from beginning of AttachList table-in Coverage Index order",
            ),
        ],
    ),
    (
        "AttachPoint",
        [
            FieldSpec(
                "uint16",
                "PointCount",
                description="Number of attachment points on this glyph",
            ),
            FieldSpec(
                "uint16",
                "PointIndex",
                repeat="PointCount",
                aux=0,
                description="Array of contour point indices -in increasing numerical order",
            ),
        ],
    ),
    (
        "LigCaretList",
        [
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table - from beginning of LigCaretList table",
            ),
            FieldSpec(
                "uint16", "LigGlyphCount", description="Number of ligature glyphs"
            ),
            FieldSpec(
                "Offset",
                "LigGlyph",
                repeat="LigGlyphCount",
                aux=0,
                description="Array of offsets to LigGlyph tables-from beginning of LigCaretList table-in Coverage Index order",
            ),
        ],
    ),
    (
        "LigGlyph",
        [
            FieldSpec(
                "uint16",
                "CaretCount",
                description="Number of CaretValues for this ligature (components - 1)",
            ),
            FieldSpec(
                "Offset",
                "CaretValue",
                repeat="CaretCount",
                aux=0,
                description="Array of offsets to CaretValue tables-from beginning of LigGlyph table-in increasing coordinate order",
            ),
        ],
    ),
    (
        "CaretValueFormat1",
        [
            FieldSpec(
                "uint16", "CaretValueFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "int16", "Coordinate", description="X or Y value, in design units"
            ),
        ],
    ),
    (
        "CaretValueFormat2",
        [
            FieldSpec(
                "uint16", "CaretValueFormat", description="Format identifier-format = 2"
            ),
            FieldSpec(
                "uint16", "CaretValuePoint", description="Contour point index on glyph"
            ),
        ],
    ),
    (
        "CaretValueFormat3",
        [
            FieldSpec(
                "uint16", "CaretValueFormat", description="Format identifier-format = 3"
            ),
            FieldSpec(
                "int16", "Coordinate", description="X or Y value, in design units"
            ),
            FieldSpec(
                "Offset",
                "DeviceTable",
                description="Offset to Device table for X or Y value-from beginning of CaretValue table",
            ),
        ],
    ),
    (
        "MarkGlyphSetsDef",
        [
            FieldSpec(
                "uint16", "MarkSetTableFormat", description="Format identifier == 1"
            ),
            FieldSpec(
                "uint16", "MarkSetCount", description="Number of mark sets defined"
            ),
            FieldSpec(
                "LOffset",
                "Coverage",
                repeat="MarkSetCount",
                aux=0,
                description="Array of offsets to mark set coverage tables.",
            ),
        ],
    ),
    #
    # base
    #
    (
        "BASE",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version of the BASE table-initially 0x00010000",
            ),
            FieldSpec(
                "Offset",
                "HorizAxis",
                description="Offset to horizontal Axis table-from beginning of BASE table-may be NULL",
            ),
            FieldSpec(
                "Offset",
                "VertAxis",
                description="Offset to vertical Axis table-from beginning of BASE table-may be NULL",
            ),
            FieldSpec(
                "LOffset",
                "VarStore",
                aux="Version >= 0x00010001",
                description="Offset to variation store (may be NULL)",
            ),
        ],
    ),
    (
        "Axis",
        [
            FieldSpec(
                "Offset",
                "BaseTagList",
                description="Offset to BaseTagList table-from beginning of Axis table-may be NULL",
            ),
            FieldSpec(
                "Offset",
                "BaseScriptList",
                description="Offset to BaseScriptList table-from beginning of Axis table",
            ),
        ],
    ),
    (
        "BaseTagList",
        [
            FieldSpec(
                "uint16",
                "BaseTagCount",
                description="Number of baseline identification tags in this text direction-may be zero (0)",
            ),
            FieldSpec(
                "Tag",
                "BaselineTag",
                repeat="BaseTagCount",
                aux=0,
                description="Array of 4-byte baseline identification tags-must be in alphabetical order",
            ),
        ],
    ),
    (
        "BaseScriptList",
        [
            FieldSpec(
                "uint16",
                "BaseScriptCount",
                description="Number of BaseScriptRecords defined",
            ),
            FieldSpec(
                "struct",
                "BaseScriptRecord",
                repeat="BaseScriptCount",
                aux=0,
                description="Array of BaseScriptRecords-in alphabetical order by BaseScriptTag",
            ),
        ],
    ),
    (
        "BaseScriptRecord",
        [
            FieldSpec(
                "Tag", "BaseScriptTag", description="4-byte script identification tag"
            ),
            FieldSpec(
                "Offset",
                "BaseScript",
                description="Offset to BaseScript table-from beginning of BaseScriptList",
            ),
        ],
    ),
    (
        "BaseScript",
        [
            FieldSpec(
                "Offset",
                "BaseValues",
                description="Offset to BaseValues table-from beginning of BaseScript table-may be NULL",
            ),
            FieldSpec(
                "Offset",
                "DefaultMinMax",
                description="Offset to MinMax table- from beginning of BaseScript table-may be NULL",
            ),
            FieldSpec(
                "uint16",
                "BaseLangSysCount",
                description="Number of BaseLangSysRecords defined-may be zero (0)",
            ),
            FieldSpec(
                "struct",
                "BaseLangSysRecord",
                repeat="BaseLangSysCount",
                aux=0,
                description="Array of BaseLangSysRecords-in alphabetical order by BaseLangSysTag",
            ),
        ],
    ),
    (
        "BaseLangSysRecord",
        [
            FieldSpec(
                "Tag",
                "BaseLangSysTag",
                description="4-byte language system identification tag",
            ),
            FieldSpec(
                "Offset",
                "MinMax",
                description="Offset to MinMax table-from beginning of BaseScript table",
            ),
        ],
    ),
    (
        "BaseValues",
        [
            FieldSpec(
                "uint16",
                "DefaultIndex",
                description="Index number of default baseline for this script-equals index position of baseline tag in BaselineArray of the BaseTagList",
            ),
            FieldSpec(
                "uint16",
                "BaseCoordCount",
                description="Number of BaseCoord tables defined-should equal BaseTagCount in the BaseTagList",
            ),
            FieldSpec(
                "Offset",
                "BaseCoord",
                repeat="BaseCoordCount",
                aux=0,
                description="Array of offsets to BaseCoord-from beginning of BaseValues table-order matches BaselineTag array in the BaseTagList",
            ),
        ],
    ),
    (
        "MinMax",
        [
            FieldSpec(
                "Offset",
                "MinCoord",
                description="Offset to BaseCoord table-defines minimum extent value-from the beginning of MinMax table-may be NULL",
            ),
            FieldSpec(
                "Offset",
                "MaxCoord",
                description="Offset to BaseCoord table-defines maximum extent value-from the beginning of MinMax table-may be NULL",
            ),
            FieldSpec(
                "uint16",
                "FeatMinMaxCount",
                description="Number of FeatMinMaxRecords-may be zero (0)",
            ),
            FieldSpec(
                "struct",
                "FeatMinMaxRecord",
                repeat="FeatMinMaxCount",
                aux=0,
                description="Array of FeatMinMaxRecords-in alphabetical order, by FeatureTableTag",
            ),
        ],
    ),
    (
        "FeatMinMaxRecord",
        [
            FieldSpec(
                "Tag",
                "FeatureTableTag",
                description="4-byte feature identification tag-must match FeatureTag in FeatureList",
            ),
            FieldSpec(
                "Offset",
                "MinCoord",
                description="Offset to BaseCoord table-defines minimum extent value-from beginning of MinMax table-may be NULL",
            ),
            FieldSpec(
                "Offset",
                "MaxCoord",
                description="Offset to BaseCoord table-defines maximum extent value-from beginning of MinMax table-may be NULL",
            ),
        ],
    ),
    (
        "BaseCoordFormat1",
        [
            FieldSpec(
                "uint16", "BaseCoordFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "int16", "Coordinate", description="X or Y value, in design units"
            ),
        ],
    ),
    (
        "BaseCoordFormat2",
        [
            FieldSpec(
                "uint16", "BaseCoordFormat", description="Format identifier-format = 2"
            ),
            FieldSpec(
                "int16", "Coordinate", description="X or Y value, in design units"
            ),
            FieldSpec(
                "GlyphID", "ReferenceGlyph", description="GlyphID of control glyph"
            ),
            FieldSpec(
                "uint16",
                "BaseCoordPoint",
                description="Index of contour point on the ReferenceGlyph",
            ),
        ],
    ),
    (
        "BaseCoordFormat3",
        [
            FieldSpec(
                "uint16", "BaseCoordFormat", description="Format identifier-format = 3"
            ),
            FieldSpec(
                "int16", "Coordinate", description="X or Y value, in design units"
            ),
            FieldSpec(
                "Offset",
                "DeviceTable",
                description="Offset to Device table for X or Y value",
            ),
        ],
    ),
    #
    # jstf
    #
    (
        "JSTF",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version of the JSTF table-initially set to 0x00010000",
            ),
            FieldSpec(
                "uint16",
                "JstfScriptCount",
                description="Number of JstfScriptRecords in this table",
            ),
            FieldSpec(
                "struct",
                "JstfScriptRecord",
                repeat="JstfScriptCount",
                aux=0,
                description="Array of JstfScriptRecords-in alphabetical order, by JstfScriptTag",
            ),
        ],
    ),
    (
        "JstfScriptRecord",
        [
            FieldSpec(
                "Tag", "JstfScriptTag", description="4-byte JstfScript identification"
            ),
            FieldSpec(
                "Offset",
                "JstfScript",
                description="Offset to JstfScript table-from beginning of JSTF Header",
            ),
        ],
    ),
    (
        "JstfScript",
        [
            FieldSpec(
                "Offset",
                "ExtenderGlyph",
                description="Offset to ExtenderGlyph table-from beginning of JstfScript table-may be NULL",
            ),
            FieldSpec(
                "Offset",
                "DefJstfLangSys",
                description="Offset to Default JstfLangSys table-from beginning of JstfScript table-may be NULL",
            ),
            FieldSpec(
                "uint16",
                "JstfLangSysCount",
                description="Number of JstfLangSysRecords in this table- may be zero (0)",
            ),
            FieldSpec(
                "struct",
                "JstfLangSysRecord",
                repeat="JstfLangSysCount",
                aux=0,
                description="Array of JstfLangSysRecords-in alphabetical order, by JstfLangSysTag",
            ),
        ],
    ),
    (
        "JstfLangSysRecord",
        [
            FieldSpec(
                "Tag", "JstfLangSysTag", description="4-byte JstfLangSys identifier"
            ),
            FieldSpec(
                "Offset",
                "JstfLangSys",
                description="Offset to JstfLangSys table-from beginning of JstfScript table",
            ),
        ],
    ),
    (
        "ExtenderGlyph",
        [
            FieldSpec(
                "uint16",
                "GlyphCount",
                description="Number of Extender Glyphs in this script",
            ),
            FieldSpec(
                "GlyphID",
                "ExtenderGlyph",
                repeat="GlyphCount",
                aux=0,
                description="GlyphIDs-in increasing numerical order",
            ),
        ],
    ),
    (
        "JstfLangSys",
        [
            FieldSpec(
                "uint16",
                "JstfPriorityCount",
                description="Number of JstfPriority tables",
            ),
            FieldSpec(
                "Offset",
                "JstfPriority",
                repeat="JstfPriorityCount",
                aux=0,
                description="Array of offsets to JstfPriority tables-from beginning of JstfLangSys table-in priority order",
            ),
        ],
    ),
    (
        "JstfPriority",
        [
            FieldSpec(
                "Offset",
                "ShrinkageEnableGSUB",
                description="Offset to Shrinkage Enable JstfGSUBModList table-from beginning of JstfPriority table-may be NULL",
            ),
            FieldSpec(
                "Offset",
                "ShrinkageDisableGSUB",
                description="Offset to Shrinkage Disable JstfGSUBModList table-from beginning of JstfPriority table-may be NULL",
            ),
            FieldSpec(
                "Offset",
                "ShrinkageEnableGPOS",
                description="Offset to Shrinkage Enable JstfGPOSModList table-from beginning of JstfPriority table-may be NULL",
            ),
            FieldSpec(
                "Offset",
                "ShrinkageDisableGPOS",
                description="Offset to Shrinkage Disable JstfGPOSModList table-from beginning of JstfPriority table-may be NULL",
            ),
            FieldSpec(
                "Offset",
                "ShrinkageJstfMax",
                description="Offset to Shrinkage JstfMax table-from beginning of JstfPriority table -may be NULL",
            ),
            FieldSpec(
                "Offset",
                "ExtensionEnableGSUB",
                description="Offset to Extension Enable JstfGSUBModList table-may be NULL",
            ),
            FieldSpec(
                "Offset",
                "ExtensionDisableGSUB",
                description="Offset to Extension Disable JstfGSUBModList table-from beginning of JstfPriority table-may be NULL",
            ),
            FieldSpec(
                "Offset",
                "ExtensionEnableGPOS",
                description="Offset to Extension Enable JstfGSUBModList table-may be NULL",
            ),
            FieldSpec(
                "Offset",
                "ExtensionDisableGPOS",
                description="Offset to Extension Disable JstfGSUBModList table-from beginning of JstfPriority table-may be NULL",
            ),
            FieldSpec(
                "Offset",
                "ExtensionJstfMax",
                description="Offset to Extension JstfMax table-from beginning of JstfPriority table -may be NULL",
            ),
        ],
    ),
    (
        "JstfGSUBModList",
        [
            FieldSpec(
                "uint16",
                "LookupCount",
                description="Number of lookups for this modification",
            ),
            FieldSpec(
                "uint16",
                "GSUBLookupIndex",
                repeat="LookupCount",
                aux=0,
                description="Array of LookupIndex identifiers in GSUB-in increasing numerical order",
            ),
        ],
    ),
    (
        "JstfGPOSModList",
        [
            FieldSpec(
                "uint16",
                "LookupCount",
                description="Number of lookups for this modification",
            ),
            FieldSpec(
                "uint16",
                "GPOSLookupIndex",
                repeat="LookupCount",
                aux=0,
                description="Array of LookupIndex identifiers in GPOS-in increasing numerical order",
            ),
        ],
    ),
    (
        "JstfMax",
        [
            FieldSpec(
                "uint16",
                "LookupCount",
                description="Number of lookup Indices for this modification",
            ),
            FieldSpec(
                "Offset",
                "Lookup",
                repeat="LookupCount",
                aux=0,
                description="Array of offsets to GPOS-type lookup tables-from beginning of JstfMax table-in design order",
            ),
        ],
    ),
    #
    # STAT
    #
    (
        "STAT",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version of the table-initially set to 0x00010000, currently 0x00010002.",
            ),
            FieldSpec(
                "uint16",
                "DesignAxisRecordSize",
                description="Size in bytes of each design axis record",
            ),
            FieldSpec(
                "uint16", "DesignAxisCount", description="Number of design axis records"
            ),
            FieldSpec(
                "LOffsetTo(AxisRecordArray)",
                "DesignAxisRecord",
                description="Offset in bytes from the beginning of the STAT table to the start of the design axes array",
            ),
            FieldSpec(
                "uint16", "AxisValueCount", description="Number of axis value tables"
            ),
            FieldSpec(
                "LOffsetTo(AxisValueArray)",
                "AxisValueArray",
                description="Offset in bytes from the beginning of the STAT table to the start of the axes value offset array",
            ),
            FieldSpec(
                "NameID",
                "ElidedFallbackNameID",
                aux="Version >= 0x00010001",
                description="NameID to use when all style attributes are elided.",
            ),
        ],
    ),
    (
        "AxisRecordArray",
        [
            FieldSpec(
                "AxisRecord",
                "Axis",
                repeat="DesignAxisCount",
                aux=0,
                description="Axis records",
            ),
        ],
    ),
    (
        "AxisRecord",
        [
            FieldSpec(
                "Tag",
                "AxisTag",
                description="A tag identifying the axis of design variation",
            ),
            FieldSpec(
                "NameID",
                "AxisNameID",
                description='The name ID for entries in the "name" table that provide a display string for this axis',
            ),
            FieldSpec(
                "uint16",
                "AxisOrdering",
                description="A value that applications can use to determine primary sorting of face names, or for ordering of descriptors when composing family or face names",
            ),
            FieldSpec(
                "uint8",
                "MoreBytes",
                repeat="DesignAxisRecordSize",
                aux=-8,
                description="Extra bytes.  Set to empty array.",
            ),
        ],
    ),
    (
        "AxisValueArray",
        [
            FieldSpec(
                "Offset",
                "AxisValue",
                repeat="AxisValueCount",
                aux=0,
                description="Axis values",
            ),
        ],
    ),
    (
        "AxisValueFormat1",
        [
            FieldSpec("uint16", "Format", description="Format, = 1"),
            FieldSpec(
                "uint16",
                "AxisIndex",
                description="Index into the axis record array identifying the axis of design variation to which the axis value record applies.",
            ),
            FieldSpec("STATFlags", "Flags", description="Flags."),
            FieldSpec("NameID", "ValueNameID"),
            FieldSpec("Fixed", "Value"),
        ],
    ),
    (
        "AxisValueFormat2",
        [
            FieldSpec("uint16", "Format", description="Format, = 2"),
            FieldSpec(
                "uint16",
                "AxisIndex",
                description="Index into the axis record array identifying the axis of design variation to which the axis value record applies.",
            ),
            FieldSpec("STATFlags", "Flags", description="Flags."),
            FieldSpec("NameID", "ValueNameID"),
            FieldSpec("Fixed", "NominalValue"),
            FieldSpec("Fixed", "RangeMinValue"),
            FieldSpec("Fixed", "RangeMaxValue"),
        ],
    ),
    (
        "AxisValueFormat3",
        [
            FieldSpec("uint16", "Format", description="Format, = 3"),
            FieldSpec(
                "uint16",
                "AxisIndex",
                description="Index into the axis record array identifying the axis of design variation to which the axis value record applies.",
            ),
            FieldSpec("STATFlags", "Flags", description="Flags."),
            FieldSpec("NameID", "ValueNameID"),
            FieldSpec("Fixed", "Value"),
            FieldSpec("Fixed", "LinkedValue"),
        ],
    ),
    (
        "AxisValueFormat4",
        [
            FieldSpec("uint16", "Format", description="Format, = 4"),
            FieldSpec(
                "uint16",
                "AxisCount",
                description="The total number of axes contributing to this axis-values combination.",
            ),
            FieldSpec("STATFlags", "Flags", description="Flags."),
            FieldSpec("NameID", "ValueNameID"),
            FieldSpec(
                "struct",
                "AxisValueRecord",
                repeat="AxisCount",
                aux=0,
                description="Array of AxisValue records that provide the combination of axis values, one for each contributing axis. ",
            ),
        ],
    ),
    (
        "AxisValueRecord",
        [
            FieldSpec(
                "uint16",
                "AxisIndex",
                description="Index into the axis record array identifying the axis of design variation to which the axis value record applies.",
            ),
            FieldSpec(
                "Fixed",
                "Value",
                description="A numeric value for this attribute value.",
            ),
        ],
    ),
    #
    # Variation fonts
    #
    # GSUB/GPOS FeatureVariations
    (
        "FeatureVariations",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version of the table-initially set to 0x00010000",
            ),
            FieldSpec(
                "uint32",
                "FeatureVariationCount",
                description="Number of records in the FeatureVariationRecord array",
            ),
            FieldSpec(
                "struct",
                "FeatureVariationRecord",
                repeat="FeatureVariationCount",
                aux=0,
                description="Array of FeatureVariationRecord",
            ),
        ],
    ),
    (
        "FeatureVariationRecord",
        [
            FieldSpec(
                "LOffset",
                "ConditionSet",
                description="Offset to a ConditionSet table, from beginning of the FeatureVariations table.",
            ),
            FieldSpec(
                "LOffset",
                "FeatureTableSubstitution",
                description="Offset to a FeatureTableSubstitution table, from beginning of the FeatureVariations table",
            ),
        ],
    ),
    (
        "ConditionList",
        [
            FieldSpec(
                "uint32",
                "ConditionCount",
                description="Number of condition tables in the ConditionTable array",
            ),
            FieldSpec(
                "LOffset",
                "ConditionTable",
                repeat="ConditionCount",
                aux=0,
                description="Array of offset to condition tables, from the beginning of the ConditionList table.",
            ),
        ],
    ),
    (
        "ConditionSet",
        [
            FieldSpec(
                "uint16",
                "ConditionCount",
                description="Number of condition tables in the ConditionTable array",
            ),
            FieldSpec(
                "LOffset",
                "ConditionTable",
                repeat="ConditionCount",
                aux=0,
                description="Array of offset to condition tables, from the beginning of the ConditionSet table.",
            ),
        ],
    ),
    (
        "ConditionTableFormat1",
        [
            FieldSpec("uint16", "Format", description="Format, = 1"),
            FieldSpec(
                "uint16",
                "AxisIndex",
                description="Index for the variation axis within the fvar table, base 0.",
            ),
            FieldSpec(
                "F2Dot14",
                "FilterRangeMinValue",
                description="Minimum normalized axis value of the font variation instances that satisfy this condition.",
            ),
            FieldSpec(
                "F2Dot14",
                "FilterRangeMaxValue",
                description="Maximum value that satisfies this condition.",
            ),
        ],
    ),
    (
        "ConditionTableFormat2",
        [
            FieldSpec("uint16", "Format", description="Format, = 2"),
            FieldSpec(
                "int16", "DefaultValue", description="Value at default instance."
            ),
            FieldSpec(
                "uint32",
                "VarIdx",
                description="Variation index to vary the value based on current designspace location.",
            ),
        ],
    ),
    (
        "ConditionTableFormat3",
        [
            FieldSpec("uint16", "Format", description="Format, = 3"),
            FieldSpec(
                "uint8",
                "ConditionCount",
                description="Index for the variation axis within the fvar table, base 0.",
            ),
            FieldSpec(
                "Offset24",
                "ConditionTable",
                repeat="ConditionCount",
                aux=0,
                description="Array of condition tables for this conjunction (AND) expression.",
            ),
        ],
    ),
    (
        "ConditionTableFormat4",
        [
            FieldSpec("uint16", "Format", description="Format, = 4"),
            FieldSpec(
                "uint8",
                "ConditionCount",
                description="Index for the variation axis within the fvar table, base 0.",
            ),
            FieldSpec(
                "Offset24",
                "ConditionTable",
                repeat="ConditionCount",
                aux=0,
                description="Array of condition tables for this disjunction (OR) expression.",
            ),
        ],
    ),
    (
        "ConditionTableFormat5",
        [
            FieldSpec("uint16", "Format", description="Format, = 5"),
            FieldSpec("Offset24", "ConditionTable", description="Condition to negate."),
        ],
    ),
    (
        "FeatureTableSubstitution",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version of the table-initially set to 0x00010000",
            ),
            FieldSpec(
                "uint16",
                "SubstitutionCount",
                description="Number of records in the FeatureVariationRecords array",
            ),
            FieldSpec(
                "FeatureTableSubstitutionRecord",
                "SubstitutionRecord",
                repeat="SubstitutionCount",
                aux=0,
                description="Array of FeatureTableSubstitutionRecord",
            ),
        ],
    ),
    (
        "FeatureTableSubstitutionRecord",
        [
            FieldSpec(
                "uint16",
                "FeatureIndex",
                description="The feature table index to match.",
            ),
            FieldSpec(
                "LOffset",
                "Feature",
                description="Offset to an alternate feature table, from start of the FeatureTableSubstitution table.",
            ),
        ],
    ),
    # VariationStore
    (
        "VarRegionAxis",
        [
            FieldSpec("F2Dot14", "StartCoord"),
            FieldSpec("F2Dot14", "PeakCoord"),
            FieldSpec("F2Dot14", "EndCoord"),
        ],
    ),
    (
        "VarRegion",
        [
            FieldSpec("struct", "VarRegionAxis", repeat="RegionAxisCount", aux=0),
        ],
    ),
    (
        "VarRegionList",
        [
            FieldSpec("uint16", "RegionAxisCount"),
            FieldSpec("uint16", "RegionCount"),
            FieldSpec("VarRegion", "Region", repeat="RegionCount", aux=0),
        ],
    ),
    (
        "VarData",
        [
            FieldSpec("uint16", "ItemCount"),
            FieldSpec("uint16", "NumShorts"),
            FieldSpec("uint16", "VarRegionCount"),
            FieldSpec("uint16", "VarRegionIndex", repeat="VarRegionCount", aux=0),
            FieldSpec("VarDataValue", "Item", repeat="ItemCount", aux=0),
        ],
    ),
    (
        "VarStore",
        [
            FieldSpec("uint16", "Format", description="Set to 1."),
            FieldSpec("LOffset", "VarRegionList"),
            FieldSpec("uint16", "VarDataCount"),
            FieldSpec("LOffset", "VarData", repeat="VarDataCount", aux=0),
        ],
    ),
    # Variation helpers
    (
        "VarIdxMap",
        [
            FieldSpec("uint16", "EntryFormat"),  # Automatically computed
            FieldSpec("uint16", "MappingCount"),  # Automatically computed
            FieldSpec(
                "VarIdxMapValue",
                "mapping",
                repeat="",
                aux=0,
                description="Array of compressed data",
            ),
        ],
    ),
    (
        "DeltaSetIndexMapFormat0",
        [
            FieldSpec(
                "uint8", "Format", description="Format of the DeltaSetIndexMap = 0"
            ),
            FieldSpec("uint8", "EntryFormat"),  # Automatically computed
            FieldSpec("uint16", "MappingCount"),  # Automatically computed
            FieldSpec(
                "VarIdxMapValue",
                "mapping",
                repeat="",
                aux=0,
                description="Array of compressed data",
            ),
        ],
    ),
    (
        "DeltaSetIndexMapFormat1",
        [
            FieldSpec(
                "uint8", "Format", description="Format of the DeltaSetIndexMap = 1"
            ),
            FieldSpec("uint8", "EntryFormat"),  # Automatically computed
            FieldSpec("uint32", "MappingCount"),  # Automatically computed
            FieldSpec(
                "VarIdxMapValue",
                "mapping",
                repeat="",
                aux=0,
                description="Array of compressed data",
            ),
        ],
    ),
    # MultiVariationStore
    (
        "SparseVarRegionAxis",
        [
            FieldSpec("uint16", "AxisIndex"),
            FieldSpec("F2Dot14", "StartCoord"),
            FieldSpec("F2Dot14", "PeakCoord"),
            FieldSpec("F2Dot14", "EndCoord"),
        ],
    ),
    (
        "SparseVarRegion",
        [
            FieldSpec("uint16", "SparseRegionCount"),
            FieldSpec(
                "struct", "SparseVarRegionAxis", repeat="SparseRegionCount", aux=0
            ),
        ],
    ),
    (
        "SparseVarRegionList",
        [
            FieldSpec("uint16", "RegionCount"),
            FieldSpec(
                "LOffsetTo(SparseVarRegion)", "Region", repeat="RegionCount", aux=0
            ),
        ],
    ),
    (
        "MultiVarData",
        [
            FieldSpec("uint8", "Format", description="Set to 1."),
            FieldSpec("uint16", "VarRegionCount"),
            FieldSpec("uint16", "VarRegionIndex", repeat="VarRegionCount", aux=0),
            FieldSpec("TupleList", "Item", repeat="", aux=0),
        ],
    ),
    (
        "MultiVarStore",
        [
            FieldSpec("uint16", "Format", description="Set to 1."),
            FieldSpec("LOffset", "SparseVarRegionList"),
            FieldSpec("uint16", "MultiVarDataCount"),
            FieldSpec("LOffset", "MultiVarData", repeat="MultiVarDataCount", aux=0),
        ],
    ),
    # VariableComposites
    (
        "VARC",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version of the HVAR table-initially = 0x00010000",
            ),
            FieldSpec("LOffset", "Coverage"),
            FieldSpec("LOffset", "MultiVarStore", description="(may be NULL)"),
            FieldSpec("LOffset", "ConditionList", description="(may be NULL)"),
            FieldSpec("LOffset", "AxisIndicesList", description="(may be NULL)"),
            FieldSpec("LOffset", "VarCompositeGlyphs"),
        ],
    ),
    (
        "AxisIndicesList",
        [
            FieldSpec("TupleList", "Item", repeat="", aux=0),
        ],
    ),
    (
        "VarCompositeGlyphs",
        [
            FieldSpec("VarCompositeGlyphList", "VarCompositeGlyph", repeat=""),
        ],
    ),
    # Glyph advance variations
    (
        "HVAR",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version of the HVAR table-initially = 0x00010000",
            ),
            FieldSpec("LOffset", "VarStore"),
            FieldSpec("LOffsetTo(VarIdxMap)", "AdvWidthMap"),
            FieldSpec("LOffsetTo(VarIdxMap)", "LsbMap"),
            FieldSpec("LOffsetTo(VarIdxMap)", "RsbMap"),
        ],
    ),
    (
        "VVAR",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version of the VVAR table-initially = 0x00010000",
            ),
            FieldSpec("LOffset", "VarStore"),
            FieldSpec("LOffsetTo(VarIdxMap)", "AdvHeightMap"),
            FieldSpec("LOffsetTo(VarIdxMap)", "TsbMap"),
            FieldSpec("LOffsetTo(VarIdxMap)", "BsbMap"),
            FieldSpec(
                "LOffsetTo(VarIdxMap)",
                "VOrgMap",
                description="Vertical origin mapping.",
            ),
        ],
    ),
    # Font-wide metrics variations
    (
        "MetricsValueRecord",
        [
            FieldSpec(
                "Tag", "ValueTag", description="4-byte font-wide measure identifier"
            ),
            FieldSpec(
                "uint32", "VarIdx", description="Combined outer-inner variation index"
            ),
            FieldSpec(
                "uint8",
                "MoreBytes",
                repeat="ValueRecordSize",
                aux=-8,
                description="Extra bytes.  Set to empty array.",
            ),
        ],
    ),
    (
        "MVAR",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version of the MVAR table-initially = 0x00010000",
            ),
            FieldSpec("uint16", "Reserved", description="Set to 0"),
            FieldSpec("uint16", "ValueRecordSize"),
            FieldSpec("uint16", "ValueRecordCount"),
            FieldSpec("Offset", "VarStore"),
            FieldSpec(
                "MetricsValueRecord", "ValueRecord", repeat="ValueRecordCount", aux=0
            ),
        ],
    ),
    #
    # math
    #
    (
        "MATH",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version of the MATH table-initially set to 0x00010000.",
            ),
            FieldSpec(
                "Offset",
                "MathConstants",
                description="Offset to MathConstants table - from the beginning of MATH table.",
            ),
            FieldSpec(
                "Offset",
                "MathGlyphInfo",
                description="Offset to MathGlyphInfo table - from the beginning of MATH table.",
            ),
            FieldSpec(
                "Offset",
                "MathVariants",
                description="Offset to MathVariants table - from the beginning of MATH table.",
            ),
        ],
    ),
    (
        "MathValueRecord",
        [
            FieldSpec(
                "int16", "Value", description="The X or Y value in design units."
            ),
            FieldSpec(
                "Offset",
                "DeviceTable",
                description="Offset to the device table - from the beginning of parent table. May be NULL. Suggested format for device table is 1.",
            ),
        ],
    ),
    (
        "MathConstants",
        [
            FieldSpec(
                "int16",
                "ScriptPercentScaleDown",
                description="Percentage of scaling down for script level 1. Suggested value: 80%.",
            ),
            FieldSpec(
                "int16",
                "ScriptScriptPercentScaleDown",
                description="Percentage of scaling down for script level 2 (ScriptScript). Suggested value: 60%.",
            ),
            FieldSpec(
                "uint16",
                "DelimitedSubFormulaMinHeight",
                description="Minimum height required for a delimited expression to be treated as a subformula. Suggested value: normal line height x1.5.",
            ),
            FieldSpec(
                "uint16",
                "DisplayOperatorMinHeight",
                description="Minimum height of n-ary operators (such as integral and summation) for formulas in display mode.",
            ),
            FieldSpec(
                "MathValueRecord",
                "MathLeading",
                description="White space to be left between math formulas to ensure proper line spacing. For example, for applications that treat line gap as a part of line ascender, formulas with ink  going above (os2.sTypoAscender + os2.sTypoLineGap - MathLeading) or with ink going below os2.sTypoDescender will result in increasing line height.",
            ),
            FieldSpec(
                "MathValueRecord", "AxisHeight", description="Axis height of the font."
            ),
            FieldSpec(
                "MathValueRecord",
                "AccentBaseHeight",
                description="Maximum (ink) height of accent base that does not require raising the accents. Suggested: x-height of the font (os2.sxHeight) plus any possible overshots.",
            ),
            FieldSpec(
                "MathValueRecord",
                "FlattenedAccentBaseHeight",
                description="Maximum (ink) height of accent base that does not require flattening the accents. Suggested: cap height of the font (os2.sCapHeight).",
            ),
            FieldSpec(
                "MathValueRecord",
                "SubscriptShiftDown",
                description="The standard shift down applied to subscript elements. Positive for moving in the downward direction. Suggested: os2.ySubscriptYOffset.",
            ),
            FieldSpec(
                "MathValueRecord",
                "SubscriptTopMax",
                description="Maximum allowed height of the (ink) top of subscripts that does not require moving subscripts further down. Suggested: 4/5 x-height.",
            ),
            FieldSpec(
                "MathValueRecord",
                "SubscriptBaselineDropMin",
                description="Minimum allowed drop of the baseline of subscripts relative to the (ink) bottom of the base. Checked for bases that are treated as a box or extended shape. Positive for subscript baseline dropped below the base bottom.",
            ),
            FieldSpec(
                "MathValueRecord",
                "SuperscriptShiftUp",
                description="Standard shift up applied to superscript elements. Suggested: os2.ySuperscriptYOffset.",
            ),
            FieldSpec(
                "MathValueRecord",
                "SuperscriptShiftUpCramped",
                description="Standard shift of superscripts relative to the base, in cramped style.",
            ),
            FieldSpec(
                "MathValueRecord",
                "SuperscriptBottomMin",
                description="Minimum allowed height of the (ink) bottom of superscripts that does not require moving subscripts further up. Suggested: 1/4 x-height.",
            ),
            FieldSpec(
                "MathValueRecord",
                "SuperscriptBaselineDropMax",
                description="Maximum allowed drop of the baseline of superscripts relative to the (ink) top of the base. Checked for bases that are treated as a box or extended shape. Positive for superscript baseline below the base top.",
            ),
            FieldSpec(
                "MathValueRecord",
                "SubSuperscriptGapMin",
                description="Minimum gap between the superscript and subscript ink. Suggested: 4x default rule thickness.",
            ),
            FieldSpec(
                "MathValueRecord",
                "SuperscriptBottomMaxWithSubscript",
                description="The maximum level to which the (ink) bottom of superscript can be pushed to increase the gap between superscript and subscript, before subscript starts being moved down. Suggested: 4/5 x-height.",
            ),
            FieldSpec(
                "MathValueRecord",
                "SpaceAfterScript",
                description="Extra white space to be added after each subscript and superscript. Suggested: 0.5pt for a 12 pt font.",
            ),
            FieldSpec(
                "MathValueRecord",
                "UpperLimitGapMin",
                description="Minimum gap between the (ink) bottom of the upper limit, and the (ink) top of the base operator.",
            ),
            FieldSpec(
                "MathValueRecord",
                "UpperLimitBaselineRiseMin",
                description="Minimum distance between baseline of upper limit and (ink) top of the base operator.",
            ),
            FieldSpec(
                "MathValueRecord",
                "LowerLimitGapMin",
                description="Minimum gap between (ink) top of the lower limit, and (ink) bottom of the base operator.",
            ),
            FieldSpec(
                "MathValueRecord",
                "LowerLimitBaselineDropMin",
                description="Minimum distance between baseline of the lower limit and (ink) bottom of the base operator.",
            ),
            FieldSpec(
                "MathValueRecord",
                "StackTopShiftUp",
                description="Standard shift up applied to the top element of a stack.",
            ),
            FieldSpec(
                "MathValueRecord",
                "StackTopDisplayStyleShiftUp",
                description="Standard shift up applied to the top element of a stack in display style.",
            ),
            FieldSpec(
                "MathValueRecord",
                "StackBottomShiftDown",
                description="Standard shift down applied to the bottom element of a stack. Positive for moving in the downward direction.",
            ),
            FieldSpec(
                "MathValueRecord",
                "StackBottomDisplayStyleShiftDown",
                description="Standard shift down applied to the bottom element of a stack in display style. Positive for moving in the downward direction.",
            ),
            FieldSpec(
                "MathValueRecord",
                "StackGapMin",
                description="Minimum gap between (ink) bottom of the top element of a stack, and the (ink) top of the bottom element. Suggested: 3x default rule thickness.",
            ),
            FieldSpec(
                "MathValueRecord",
                "StackDisplayStyleGapMin",
                description="Minimum gap between (ink) bottom of the top element of a stack, and the (ink) top of the bottom element in display style. Suggested: 7x default rule thickness.",
            ),
            FieldSpec(
                "MathValueRecord",
                "StretchStackTopShiftUp",
                description="Standard shift up applied to the top element of the stretch stack.",
            ),
            FieldSpec(
                "MathValueRecord",
                "StretchStackBottomShiftDown",
                description="Standard shift down applied to the bottom element of the stretch stack. Positive for moving in the downward direction.",
            ),
            FieldSpec(
                "MathValueRecord",
                "StretchStackGapAboveMin",
                description="Minimum gap between the ink of the stretched element, and the (ink) bottom of the element above. Suggested: UpperLimitGapMin",
            ),
            FieldSpec(
                "MathValueRecord",
                "StretchStackGapBelowMin",
                description="Minimum gap between the ink of the stretched element, and the (ink) top of the element below. Suggested: LowerLimitGapMin.",
            ),
            FieldSpec(
                "MathValueRecord",
                "FractionNumeratorShiftUp",
                description="Standard shift up applied to the numerator.",
            ),
            FieldSpec(
                "MathValueRecord",
                "FractionNumeratorDisplayStyleShiftUp",
                description="Standard shift up applied to the numerator in display style. Suggested: StackTopDisplayStyleShiftUp.",
            ),
            FieldSpec(
                "MathValueRecord",
                "FractionDenominatorShiftDown",
                description="Standard shift down applied to the denominator. Positive for moving in the downward direction.",
            ),
            FieldSpec(
                "MathValueRecord",
                "FractionDenominatorDisplayStyleShiftDown",
                description="Standard shift down applied to the denominator in display style. Positive for moving in the downward direction. Suggested: StackBottomDisplayStyleShiftDown.",
            ),
            FieldSpec(
                "MathValueRecord",
                "FractionNumeratorGapMin",
                description="Minimum tolerated gap between the (ink) bottom of the numerator and the ink of the fraction bar. Suggested: default rule thickness",
            ),
            FieldSpec(
                "MathValueRecord",
                "FractionNumDisplayStyleGapMin",
                description="Minimum tolerated gap between the (ink) bottom of the numerator and the ink of the fraction bar in display style. Suggested: 3x default rule thickness.",
            ),
            FieldSpec(
                "MathValueRecord",
                "FractionRuleThickness",
                description="Thickness of the fraction bar. Suggested: default rule thickness.",
            ),
            FieldSpec(
                "MathValueRecord",
                "FractionDenominatorGapMin",
                description="Minimum tolerated gap between the (ink) top of the denominator and the ink of the fraction bar. Suggested: default rule thickness",
            ),
            FieldSpec(
                "MathValueRecord",
                "FractionDenomDisplayStyleGapMin",
                description="Minimum tolerated gap between the (ink) top of the denominator and the ink of the fraction bar in display style. Suggested: 3x default rule thickness.",
            ),
            FieldSpec(
                "MathValueRecord",
                "SkewedFractionHorizontalGap",
                description="Horizontal distance between the top and bottom elements of a skewed fraction.",
            ),
            FieldSpec(
                "MathValueRecord",
                "SkewedFractionVerticalGap",
                description="Vertical distance between the ink of the top and bottom elements of a skewed fraction.",
            ),
            FieldSpec(
                "MathValueRecord",
                "OverbarVerticalGap",
                description="Distance between the overbar and the (ink) top of he base. Suggested: 3x default rule thickness.",
            ),
            FieldSpec(
                "MathValueRecord",
                "OverbarRuleThickness",
                description="Thickness of overbar. Suggested: default rule thickness.",
            ),
            FieldSpec(
                "MathValueRecord",
                "OverbarExtraAscender",
                description="Extra white space reserved above the overbar. Suggested: default rule thickness.",
            ),
            FieldSpec(
                "MathValueRecord",
                "UnderbarVerticalGap",
                description="Distance between underbar and (ink) bottom of the base. Suggested: 3x default rule thickness.",
            ),
            FieldSpec(
                "MathValueRecord",
                "UnderbarRuleThickness",
                description="Thickness of underbar. Suggested: default rule thickness.",
            ),
            FieldSpec(
                "MathValueRecord",
                "UnderbarExtraDescender",
                description="Extra white space reserved below the underbar. Always positive. Suggested: default rule thickness.",
            ),
            FieldSpec(
                "MathValueRecord",
                "RadicalVerticalGap",
                description="Space between the (ink) top of the expression and the bar over it. Suggested: 1 1/4 default rule thickness.",
            ),
            FieldSpec(
                "MathValueRecord",
                "RadicalDisplayStyleVerticalGap",
                description="Space between the (ink) top of the expression and the bar over it. Suggested: default rule thickness + 1/4 x-height.",
            ),
            FieldSpec(
                "MathValueRecord",
                "RadicalRuleThickness",
                description="Thickness of the radical rule. This is the thickness of the rule in designed or constructed radical signs. Suggested: default rule thickness.",
            ),
            FieldSpec(
                "MathValueRecord",
                "RadicalExtraAscender",
                description="Extra white space reserved above the radical. Suggested: RadicalRuleThickness.",
            ),
            FieldSpec(
                "MathValueRecord",
                "RadicalKernBeforeDegree",
                description="Extra horizontal kern before the degree of a radical, if such is present. Suggested: 5/18 of em.",
            ),
            FieldSpec(
                "MathValueRecord",
                "RadicalKernAfterDegree",
                description="Negative kern after the degree of a radical, if such is present. Suggested: 10/18 of em.",
            ),
            FieldSpec(
                "uint16",
                "RadicalDegreeBottomRaisePercent",
                description="Height of the bottom of the radical degree, if such is present, in proportion to the ascender of the radical sign. Suggested: 60%.",
            ),
        ],
    ),
    (
        "MathGlyphInfo",
        [
            FieldSpec(
                "Offset",
                "MathItalicsCorrectionInfo",
                description="Offset to MathItalicsCorrectionInfo table - from the beginning of MathGlyphInfo table.",
            ),
            FieldSpec(
                "Offset",
                "MathTopAccentAttachment",
                description="Offset to MathTopAccentAttachment table - from the beginning of MathGlyphInfo table.",
            ),
            FieldSpec(
                "Offset",
                "ExtendedShapeCoverage",
                description="Offset to coverage table for Extended Shape glyphs - from the  beginning of MathGlyphInfo table. When the left or right glyph of a box is an extended shape variant, the (ink) box (and not the default position defined by values in MathConstants table) should be used for vertical positioning purposes. May be NULL.",
            ),
            FieldSpec(
                "Offset",
                "MathKernInfo",
                description="Offset to MathKernInfo table - from the beginning of MathGlyphInfo table.",
            ),
        ],
    ),
    (
        "MathItalicsCorrectionInfo",
        [
            FieldSpec(
                "Offset",
                "Coverage",
                description="Offset to Coverage table - from the beginning of MathItalicsCorrectionInfo table.",
            ),
            FieldSpec(
                "uint16",
                "ItalicsCorrectionCount",
                description="Number of italics correction values. Should coincide with the number of covered glyphs.",
            ),
            FieldSpec(
                "MathValueRecord",
                "ItalicsCorrection",
                repeat="ItalicsCorrectionCount",
                aux=0,
                description="Array of MathValueRecords defining italics correction values for each covered glyph.",
            ),
        ],
    ),
    (
        "MathTopAccentAttachment",
        [
            FieldSpec(
                "Offset",
                "TopAccentCoverage",
                description="Offset to Coverage table - from the beginning of  MathTopAccentAttachment table.",
            ),
            FieldSpec(
                "uint16",
                "TopAccentAttachmentCount",
                description="Number of top accent attachment point values. Should coincide with the number of covered glyphs",
            ),
            FieldSpec(
                "MathValueRecord",
                "TopAccentAttachment",
                repeat="TopAccentAttachmentCount",
                aux=0,
                description="Array of MathValueRecords defining top accent attachment points for each covered glyph",
            ),
        ],
    ),
    (
        "MathKernInfo",
        [
            FieldSpec(
                "Offset",
                "MathKernCoverage",
                description="Offset to Coverage table - from the beginning of the MathKernInfo table.",
            ),
            FieldSpec(
                "uint16", "MathKernCount", description="Number of MathKernInfoRecords."
            ),
            FieldSpec(
                "MathKernInfoRecord",
                "MathKernInfoRecords",
                repeat="MathKernCount",
                aux=0,
                description="Array of MathKernInfoRecords, per-glyph information for mathematical positioning of subscripts and superscripts.",
            ),
        ],
    ),
    (
        "MathKernInfoRecord",
        [
            FieldSpec(
                "Offset",
                "TopRightMathKern",
                description="Offset to MathKern table for top right corner - from the beginning of MathKernInfo table. May be NULL.",
            ),
            FieldSpec(
                "Offset",
                "TopLeftMathKern",
                description="Offset to MathKern table for the top left corner - from the beginning of MathKernInfo table. May be NULL.",
            ),
            FieldSpec(
                "Offset",
                "BottomRightMathKern",
                description="Offset to MathKern table for bottom right corner - from the beginning of MathKernInfo table. May be NULL.",
            ),
            FieldSpec(
                "Offset",
                "BottomLeftMathKern",
                description="Offset to MathKern table for bottom left corner - from the beginning of MathKernInfo table. May be NULL.",
            ),
        ],
    ),
    (
        "MathKern",
        [
            FieldSpec(
                "uint16",
                "HeightCount",
                description="Number of heights on which the kern value changes.",
            ),
            FieldSpec(
                "MathValueRecord",
                "CorrectionHeight",
                repeat="HeightCount",
                aux=0,
                description="Array of correction heights at which the kern value changes. Sorted by the height value in design units.",
            ),
            FieldSpec(
                "MathValueRecord",
                "KernValue",
                repeat="HeightCount",
                aux=1,
                description="Array of kern values corresponding to heights. First value is the kern value for all heights less or equal than the first height in this table.Last value is the value to be applied for all heights greater than the last height in this table. Negative values are interpreted as move glyphs closer to each other.",
            ),
        ],
    ),
    (
        "MathVariants",
        [
            FieldSpec(
                "uint16",
                "MinConnectorOverlap",
                description="Minimum overlap of connecting glyphs during glyph construction,  in design units.",
            ),
            FieldSpec(
                "Offset",
                "VertGlyphCoverage",
                description="Offset to Coverage table - from the beginning of MathVariants table.",
            ),
            FieldSpec(
                "Offset",
                "HorizGlyphCoverage",
                description="Offset to Coverage table - from the beginning of MathVariants table.",
            ),
            FieldSpec(
                "uint16",
                "VertGlyphCount",
                description="Number of glyphs for which information is provided for vertically growing variants.",
            ),
            FieldSpec(
                "uint16",
                "HorizGlyphCount",
                description="Number of glyphs for which information is provided for horizontally growing variants.",
            ),
            FieldSpec(
                "Offset",
                "VertGlyphConstruction",
                repeat="VertGlyphCount",
                aux=0,
                description="Array of offsets to MathGlyphConstruction tables - from the beginning of the MathVariants table, for shapes growing in vertical direction.",
            ),
            FieldSpec(
                "Offset",
                "HorizGlyphConstruction",
                repeat="HorizGlyphCount",
                aux=0,
                description="Array of offsets to MathGlyphConstruction tables - from the beginning of the MathVariants table, for shapes growing in horizontal direction.",
            ),
        ],
    ),
    (
        "MathGlyphConstruction",
        [
            FieldSpec(
                "Offset",
                "GlyphAssembly",
                description="Offset to GlyphAssembly table for this shape - from the beginning of MathGlyphConstruction table. May be NULL",
            ),
            FieldSpec(
                "uint16",
                "VariantCount",
                description="Count of glyph growing variants for this glyph.",
            ),
            FieldSpec(
                "MathGlyphVariantRecord",
                "MathGlyphVariantRecord",
                repeat="VariantCount",
                aux=0,
                description="MathGlyphVariantRecords for alternative variants of the glyphs.",
            ),
        ],
    ),
    (
        "MathGlyphVariantRecord",
        [
            FieldSpec(
                "GlyphID", "VariantGlyph", description="Glyph ID for the variant."
            ),
            FieldSpec(
                "uint16",
                "AdvanceMeasurement",
                description="Advance width/height, in design units, of the variant, in the direction of requested glyph extension.",
            ),
        ],
    ),
    (
        "GlyphAssembly",
        [
            FieldSpec(
                "MathValueRecord",
                "ItalicsCorrection",
                description="Italics correction of this GlyphAssembly. Should not depend on the assembly size.",
            ),
            FieldSpec(
                "uint16", "PartCount", description="Number of parts in this assembly."
            ),
            FieldSpec(
                "GlyphPartRecord",
                "PartRecords",
                repeat="PartCount",
                aux=0,
                description="Array of part records, from left to right and bottom to top.",
            ),
        ],
    ),
    (
        "GlyphPartRecord",
        [
            FieldSpec("GlyphID", "glyph", description="Glyph ID for the part."),
            FieldSpec(
                "uint16",
                "StartConnectorLength",
                description="Advance width/ height of the straight bar connector material, in design units, is at the beginning of the glyph, in the direction of the extension.",
            ),
            FieldSpec(
                "uint16",
                "EndConnectorLength",
                description="Advance width/ height of the straight bar connector material, in design units, is at the end of the glyph, in the direction of the extension.",
            ),
            FieldSpec(
                "uint16",
                "FullAdvance",
                description="Full advance width/height for this part, in the direction of the extension. In design units.",
            ),
            FieldSpec(
                "uint16",
                "PartFlags",
                description="Part qualifiers. PartFlags enumeration currently uses only one bit: 0x0001 fExtender: If set, the part can be skipped or repeated. 0xFFFE Reserved",
            ),
        ],
    ),
    ##
    ## Apple Advanced Typography (AAT) tables
    ##
    (
        "AATLookupSegment",
        [
            FieldSpec(
                "uint16", "lastGlyph", description="Last glyph index in this segment."
            ),
            FieldSpec(
                "uint16", "firstGlyph", description="First glyph index in this segment."
            ),
            FieldSpec(
                "uint16",
                "value",
                description="A 16-bit offset from the start of the table to the data.",
            ),
        ],
    ),
    #
    # ankr
    #
    (
        "ankr",
        [
            FieldSpec("struct", "AnchorPoints", description="Anchor points table."),
        ],
    ),
    (
        "AnchorPointsFormat0",
        [
            FieldSpec(
                "uint16",
                "Format",
                description="Format of the anchor points table, = 0.",
            ),
            FieldSpec(
                "uint16", "Flags", description="Flags. Currenty unused, set to zero."
            ),
            FieldSpec(
                "AATLookupWithDataOffset(AnchorGlyphData)",
                "Anchors",
                description="Table of with anchor overrides for each glyph.",
            ),
        ],
    ),
    (
        "AnchorGlyphData",
        [
            FieldSpec(
                "uint32",
                "AnchorPointCount",
                description="Number of anchor points for this glyph.",
            ),
            FieldSpec(
                "struct",
                "AnchorPoint",
                repeat="AnchorPointCount",
                aux=0,
                description="Individual anchor points.",
            ),
        ],
    ),
    (
        "AnchorPoint",
        [
            FieldSpec(
                "int16", "XCoordinate", description="X coordinate of this anchor point."
            ),
            FieldSpec(
                "int16", "YCoordinate", description="Y coordinate of this anchor point."
            ),
        ],
    ),
    #
    # bsln
    #
    (
        "bsln",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version number of the AAT baseline table (0x00010000 for the initial version).",
            ),
            FieldSpec("struct", "Baseline", description="Baseline table."),
        ],
    ),
    (
        "BaselineFormat0",
        [
            FieldSpec(
                "uint16", "Format", description="Format of the baseline table, = 0."
            ),
            FieldSpec(
                "uint16",
                "DefaultBaseline",
                description="Default baseline value for all glyphs. This value can be from 0 through 31.",
            ),
            FieldSpec(
                "uint16",
                "Delta",
                repeat=32,
                aux=0,
                description="These are the FUnit distance deltas from the font’s natural baseline to the other baselines used in the font. A total of 32 deltas must be assigned.",
            ),
        ],
    ),
    (
        "BaselineFormat1",
        [
            FieldSpec(
                "uint16", "Format", description="Format of the baseline table, = 1."
            ),
            FieldSpec(
                "uint16",
                "DefaultBaseline",
                description="Default baseline value for all glyphs. This value can be from 0 through 31.",
            ),
            FieldSpec(
                "uint16",
                "Delta",
                repeat=32,
                aux=0,
                description="These are the FUnit distance deltas from the font’s natural baseline to the other baselines used in the font. A total of 32 deltas must be assigned.",
            ),
            FieldSpec(
                "AATLookup(uint16)",
                "BaselineValues",
                description="Lookup table that maps glyphs to their baseline values.",
            ),
        ],
    ),
    (
        "BaselineFormat2",
        [
            FieldSpec(
                "uint16", "Format", description="Format of the baseline table, = 1."
            ),
            FieldSpec(
                "uint16",
                "DefaultBaseline",
                description="Default baseline value for all glyphs. This value can be from 0 through 31.",
            ),
            FieldSpec(
                "GlyphID",
                "StandardGlyph",
                description="Glyph index of the glyph in this font to be used to set the baseline values. This glyph must contain a set of control points (whose numbers are contained in the following field) that determines baseline distances.",
            ),
            FieldSpec(
                "uint16",
                "ControlPoint",
                repeat=32,
                aux=0,
                description="Array of 32 control point numbers, associated with the standard glyph. A value of 0xFFFF means there is no corresponding control point in the standard glyph.",
            ),
        ],
    ),
    (
        "BaselineFormat3",
        [
            FieldSpec(
                "uint16", "Format", description="Format of the baseline table, = 1."
            ),
            FieldSpec(
                "uint16",
                "DefaultBaseline",
                description="Default baseline value for all glyphs. This value can be from 0 through 31.",
            ),
            FieldSpec(
                "GlyphID",
                "StandardGlyph",
                description="Glyph index of the glyph in this font to be used to set the baseline values. This glyph must contain a set of control points (whose numbers are contained in the following field) that determines baseline distances.",
            ),
            FieldSpec(
                "uint16",
                "ControlPoint",
                repeat=32,
                aux=0,
                description="Array of 32 control point numbers, associated with the standard glyph. A value of 0xFFFF means there is no corresponding control point in the standard glyph.",
            ),
            FieldSpec(
                "AATLookup(uint16)",
                "BaselineValues",
                description="Lookup table that maps glyphs to their baseline values.",
            ),
        ],
    ),
    #
    # cidg
    #
    (
        "cidg",
        [
            FieldSpec(
                "struct", "CIDGlyphMapping", description="CID-to-glyph mapping table."
            ),
        ],
    ),
    (
        "CIDGlyphMappingFormat0",
        [
            FieldSpec(
                "uint16",
                "Format",
                description="Format of the CID-to-glyph mapping table, = 0.",
            ),
            FieldSpec(
                "uint16", "DataFormat", description="Currenty unused, set to zero."
            ),
            FieldSpec(
                "uint32", "StructLength", description="Size of the table in bytes."
            ),
            FieldSpec("uint16", "Registry", description="The registry ID."),
            FieldSpec(
                "char64",
                "RegistryName",
                description="The registry name in ASCII; unused bytes should be set to 0.",
            ),
            FieldSpec("uint16", "Order", description="The order ID."),
            FieldSpec(
                "char64",
                "OrderName",
                description="The order name in ASCII; unused bytes should be set to 0.",
            ),
            FieldSpec(
                "uint16", "SupplementVersion", description="The supplement version."
            ),
            FieldSpec(
                "CIDGlyphMap",
                "Mapping",
                description="A mapping from CIDs to the glyphs in the font, starting with CID 0. If a CID from the identified collection has no glyph in the font, 0xFFFF is used",
            ),
        ],
    ),
    #
    # feat
    #
    (
        "feat",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version of the feat table-initially set to 0x00010000.",
            ),
            FieldSpec("FeatureNames", "FeatureNames", description="The feature names."),
        ],
    ),
    (
        "FeatureNames",
        [
            FieldSpec(
                "uint16",
                "FeatureNameCount",
                description="Number of entries in the feature name array.",
            ),
            FieldSpec("uint16", "Reserved1", description="Reserved (set to zero)."),
            FieldSpec("uint32", "Reserved2", description="Reserved (set to zero)."),
            FieldSpec(
                "FeatureName",
                "FeatureName",
                repeat="FeatureNameCount",
                aux=0,
                description="The feature name array.",
            ),
        ],
    ),
    (
        "FeatureName",
        [
            FieldSpec("uint16", "FeatureType", description="Feature type."),
            FieldSpec(
                "uint16",
                "SettingsCount",
                description="The number of records in the setting name array.",
            ),
            FieldSpec(
                "LOffset",
                "Settings",
                description="Offset to setting table for this feature.",
            ),
            FieldSpec(
                "uint16",
                "FeatureFlags",
                description="Single-bit flags associated with the feature type.",
            ),
            FieldSpec(
                "NameID",
                "FeatureNameID",
                description="The name table index for the feature name.",
            ),
        ],
    ),
    (
        "Settings",
        [
            FieldSpec(
                "Setting",
                "Setting",
                repeat="SettingsCount",
                aux=0,
                description="The setting array.",
            ),
        ],
    ),
    (
        "Setting",
        [
            FieldSpec("uint16", "SettingValue", description="The setting."),
            FieldSpec(
                "NameID",
                "SettingNameID",
                description="The name table index for the setting name.",
            ),
        ],
    ),
    #
    # gcid
    #
    (
        "gcid",
        [
            FieldSpec(
                "struct", "GlyphCIDMapping", description="Glyph to CID mapping table."
            ),
        ],
    ),
    (
        "GlyphCIDMappingFormat0",
        [
            FieldSpec(
                "uint16",
                "Format",
                description="Format of the glyph-to-CID mapping table, = 0.",
            ),
            FieldSpec(
                "uint16", "DataFormat", description="Currenty unused, set to zero."
            ),
            FieldSpec(
                "uint32", "StructLength", description="Size of the table in bytes."
            ),
            FieldSpec("uint16", "Registry", description="The registry ID."),
            FieldSpec(
                "char64",
                "RegistryName",
                description="The registry name in ASCII; unused bytes should be set to 0.",
            ),
            FieldSpec("uint16", "Order", description="The order ID."),
            FieldSpec(
                "char64",
                "OrderName",
                description="The order name in ASCII; unused bytes should be set to 0.",
            ),
            FieldSpec(
                "uint16", "SupplementVersion", description="The supplement version."
            ),
            FieldSpec(
                "GlyphCIDMap",
                "Mapping",
                description="The CIDs for the glyphs in the font, starting with glyph 0. If a glyph does not correspond to a CID in the identified collection, 0xFFFF is used",
            ),
        ],
    ),
    #
    # lcar
    #
    (
        "lcar",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version number of the ligature caret table (0x00010000 for the initial version).",
            ),
            FieldSpec("struct", "LigatureCarets", description="Ligature carets table."),
        ],
    ),
    (
        "LigatureCaretsFormat0",
        [
            FieldSpec(
                "uint16",
                "Format",
                description="Format of the ligature caret table. Format 0 indicates division points are distances in font units, Format 1 indicates division points are indexes of control points.",
            ),
            FieldSpec(
                "AATLookup(LigCaretDistances)",
                "Carets",
                description="Lookup table associating ligature glyphs with their caret positions, in font unit distances.",
            ),
        ],
    ),
    (
        "LigatureCaretsFormat1",
        [
            FieldSpec(
                "uint16",
                "Format",
                description="Format of the ligature caret table. Format 0 indicates division points are distances in font units, Format 1 indicates division points are indexes of control points.",
            ),
            FieldSpec(
                "AATLookup(LigCaretPoints)",
                "Carets",
                description="Lookup table associating ligature glyphs with their caret positions, as control points.",
            ),
        ],
    ),
    (
        "LigCaretDistances",
        [
            FieldSpec(
                "uint16", "DivsionPointCount", description="Number of division points."
            ),
            FieldSpec(
                "int16",
                "DivisionPoint",
                repeat="DivsionPointCount",
                aux=0,
                description="Distance in font units through which a subdivision is made orthogonally to the baseline.",
            ),
        ],
    ),
    (
        "LigCaretPoints",
        [
            FieldSpec(
                "uint16", "DivsionPointCount", description="Number of division points."
            ),
            FieldSpec(
                "int16",
                "DivisionPoint",
                repeat="DivsionPointCount",
                aux=0,
                description="The number of the control point through which a subdivision is made orthogonally to the baseline.",
            ),
        ],
    ),
    #
    # mort
    #
    (
        "mort",
        [
            FieldSpec("Version", "Version", description="Version of the mort table."),
            FieldSpec(
                "uint32",
                "MorphChainCount",
                description="Number of metamorphosis chains.",
            ),
            FieldSpec(
                "MortChain",
                "MorphChain",
                repeat="MorphChainCount",
                aux=0,
                description="Array of metamorphosis chains.",
            ),
        ],
    ),
    (
        "MortChain",
        [
            FieldSpec(
                "Flags32",
                "DefaultFlags",
                description="The default specification for subtables.",
            ),
            FieldSpec(
                "uint32",
                "StructLength",
                description="Total byte count, including this header; must be a multiple of 4.",
            ),
            FieldSpec(
                "uint16",
                "MorphFeatureCount",
                description="Number of metamorphosis feature entries.",
            ),
            FieldSpec(
                "uint16",
                "MorphSubtableCount",
                description="The number of subtables in the chain.",
            ),
            FieldSpec(
                "struct",
                "MorphFeature",
                repeat="MorphFeatureCount",
                aux=0,
                description="Array of metamorphosis features.",
            ),
            FieldSpec(
                "MortSubtable",
                "MorphSubtable",
                repeat="MorphSubtableCount",
                aux=0,
                description="Array of metamorphosis subtables.",
            ),
        ],
    ),
    (
        "MortSubtable",
        [
            FieldSpec(
                "uint16",
                "StructLength",
                description="Total subtable length, including this header.",
            ),
            FieldSpec(
                "uint8",
                "CoverageFlags",
                description="Most significant byte of coverage flags.",
            ),
            FieldSpec("uint8", "MorphType", description="Subtable type."),
            FieldSpec(
                "Flags32",
                "SubFeatureFlags",
                description="The 32-bit mask identifying which subtable this is (the subtable being executed if the AND of this value and the processed defaultFlags is nonzero).",
            ),
            FieldSpec("SubStruct", "SubStruct", description="SubTable."),
        ],
    ),
    #
    # morx
    #
    (
        "morx",
        [
            FieldSpec("uint16", "Version", description="Version of the morx table."),
            FieldSpec("uint16", "Reserved", description="Reserved (set to zero)."),
            FieldSpec(
                "uint32",
                "MorphChainCount",
                description="Number of extended metamorphosis chains.",
            ),
            FieldSpec(
                "MorxChain",
                "MorphChain",
                repeat="MorphChainCount",
                aux=0,
                description="Array of extended metamorphosis chains.",
            ),
        ],
    ),
    (
        "MorxChain",
        [
            FieldSpec(
                "Flags32",
                "DefaultFlags",
                description="The default specification for subtables.",
            ),
            FieldSpec(
                "uint32",
                "StructLength",
                description="Total byte count, including this header; must be a multiple of 4.",
            ),
            FieldSpec(
                "uint32",
                "MorphFeatureCount",
                description="Number of feature subtable entries.",
            ),
            FieldSpec(
                "uint32",
                "MorphSubtableCount",
                description="The number of subtables in the chain.",
            ),
            FieldSpec(
                "MorphFeature",
                "MorphFeature",
                repeat="MorphFeatureCount",
                aux=0,
                description="Array of metamorphosis features.",
            ),
            FieldSpec(
                "MorxSubtable",
                "MorphSubtable",
                repeat="MorphSubtableCount",
                aux=0,
                description="Array of extended metamorphosis subtables.",
            ),
        ],
    ),
    (
        "MorphFeature",
        [
            FieldSpec("uint16", "FeatureType", description="The type of feature."),
            FieldSpec(
                "uint16",
                "FeatureSetting",
                description="The feature's setting (aka selector).",
            ),
            FieldSpec(
                "Flags32",
                "EnableFlags",
                description="Flags for the settings that this feature and setting enables.",
            ),
            FieldSpec(
                "Flags32",
                "DisableFlags",
                description="Complement of flags for the settings that this feature and setting disable.",
            ),
        ],
    ),
    # Apple TrueType Reference Manual, chapter “The ‘morx’ table”,
    # section “Metamorphosis Subtables”.
    # https://developer.apple.com/fonts/TrueType-Reference-Manual/RM06/Chap6morx.html
    (
        "MorxSubtable",
        [
            FieldSpec(
                "uint32",
                "StructLength",
                description="Total subtable length, including this header.",
            ),
            FieldSpec(
                "uint8",
                "CoverageFlags",
                description="Most significant byte of coverage flags.",
            ),
            FieldSpec("uint16", "Reserved", description="Unused."),
            FieldSpec("uint8", "MorphType", description="Subtable type."),
            FieldSpec(
                "Flags32",
                "SubFeatureFlags",
                description="The 32-bit mask identifying which subtable this is (the subtable being executed if the AND of this value and the processed defaultFlags is nonzero).",
            ),
            FieldSpec("SubStruct", "SubStruct", description="SubTable."),
        ],
    ),
    (
        "StateHeader",
        [
            FieldSpec(
                "uint32",
                "ClassCount",
                description="Number of classes, which is the number of 16-bit entry indices in a single line in the state array.",
            ),
            FieldSpec(
                "uint32",
                "MorphClass",
                description="Offset from the start of this state table header to the start of the class table.",
            ),
            FieldSpec(
                "uint32",
                "StateArrayOffset",
                description="Offset from the start of this state table header to the start of the state array.",
            ),
            FieldSpec(
                "uint32",
                "EntryTableOffset",
                description="Offset from the start of this state table header to the start of the entry table.",
            ),
        ],
    ),
    (
        "RearrangementMorph",
        [
            FieldSpec(
                "STXHeader(RearrangementMorphAction)",
                "StateTable",
                description="Finite-state transducer table for indic rearrangement.",
            ),
        ],
    ),
    (
        "ContextualMorph",
        [
            FieldSpec(
                "STXHeader(ContextualMorphAction)",
                "StateTable",
                description="Finite-state transducer for contextual glyph substitution.",
            ),
        ],
    ),
    (
        "LigatureMorph",
        [
            FieldSpec(
                "STXHeader(LigatureMorphAction)",
                "StateTable",
                description="Finite-state transducer for ligature substitution.",
            ),
        ],
    ),
    (
        "NoncontextualMorph",
        [
            FieldSpec(
                "AATLookup(GlyphID)",
                "Substitution",
                description="The noncontextual glyph substitution table.",
            ),
        ],
    ),
    (
        "InsertionMorph",
        [
            FieldSpec(
                "STXHeader(InsertionMorphAction)",
                "StateTable",
                description="Finite-state transducer for glyph insertion.",
            ),
        ],
    ),
    (
        "MorphClass",
        [
            FieldSpec(
                "uint16",
                "FirstGlyph",
                description="Glyph index of the first glyph in the class table.",
            ),
            # ('uint16', 'GlyphCount', None, None, 'Number of glyphs in class table.'),
            # ('uint8', 'GlyphClass', 'GlyphCount', 0, 'The class codes (indexed by glyph index minus firstGlyph). Class codes range from 0 to the value of stateSize minus 1.'),
        ],
    ),
    # If the 'morx' table version is 3 or greater, then the last subtable in the chain is followed by a subtableGlyphCoverageArray, as described below.
    # 		('Offset', 'MarkGlyphSetsDef', None, 'round(Version*0x10000) >= 0x00010002', 'Offset to the table of mark set definitions-from beginning of GDEF header (may be NULL)'),
    #
    # prop
    #
    (
        "prop",
        [
            FieldSpec(
                "Fixed",
                "Version",
                description="Version number of the AAT glyphs property table. Version 1.0 is the initial table version. Version 2.0, which is recognized by macOS 8.5 and later, adds support for the “attaches on right” bit. Version 3.0, which gets recognized by macOS X and iOS, adds support for the additional directional properties defined in Unicode 3.0.",
            ),
            FieldSpec("struct", "GlyphProperties", description="Glyph properties."),
        ],
    ),
    (
        "GlyphPropertiesFormat0",
        [
            FieldSpec("uint16", "Format", description="Format, = 0."),
            FieldSpec(
                "uint16",
                "DefaultProperties",
                description="Default properties applied to a glyph. Since there is no lookup table in prop format 0, the default properties get applied to every glyph in the font.",
            ),
        ],
    ),
    (
        "GlyphPropertiesFormat1",
        [
            FieldSpec("uint16", "Format", description="Format, = 1."),
            FieldSpec(
                "uint16",
                "DefaultProperties",
                description="Default properties applied to a glyph if that glyph is not present in the Properties lookup table.",
            ),
            FieldSpec(
                "AATLookup(uint16)",
                "Properties",
                description="Lookup data associating glyphs with their properties.",
            ),
        ],
    ),
    #
    # opbd
    #
    (
        "opbd",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version number of the optical bounds table (0x00010000 for the initial version).",
            ),
            FieldSpec("struct", "OpticalBounds", description="Optical bounds table."),
        ],
    ),
    (
        "OpticalBoundsFormat0",
        [
            FieldSpec(
                "uint16",
                "Format",
                description="Format of the optical bounds table, = 0.",
            ),
            FieldSpec(
                "AATLookup(OpticalBoundsDeltas)",
                "OpticalBoundsDeltas",
                description="Lookup table associating glyphs with their optical bounds, given as deltas in font units.",
            ),
        ],
    ),
    (
        "OpticalBoundsFormat1",
        [
            FieldSpec(
                "uint16",
                "Format",
                description="Format of the optical bounds table, = 1.",
            ),
            FieldSpec(
                "AATLookup(OpticalBoundsPoints)",
                "OpticalBoundsPoints",
                description="Lookup table associating glyphs with their optical bounds, given as references to control points.",
            ),
        ],
    ),
    (
        "OpticalBoundsDeltas",
        [
            FieldSpec(
                "int16",
                "Left",
                description="Delta value for the left-side optical edge.",
            ),
            FieldSpec(
                "int16", "Top", description="Delta value for the top-side optical edge."
            ),
            FieldSpec(
                "int16",
                "Right",
                description="Delta value for the right-side optical edge.",
            ),
            FieldSpec(
                "int16",
                "Bottom",
                description="Delta value for the bottom-side optical edge.",
            ),
        ],
    ),
    (
        "OpticalBoundsPoints",
        [
            FieldSpec(
                "int16",
                "Left",
                description="Control point index for the left-side optical edge, or -1 if this glyph has none.",
            ),
            FieldSpec(
                "int16",
                "Top",
                description="Control point index for the top-side optical edge, or -1 if this glyph has none.",
            ),
            FieldSpec(
                "int16",
                "Right",
                description="Control point index for the right-side optical edge, or -1 if this glyph has none.",
            ),
            FieldSpec(
                "int16",
                "Bottom",
                description="Control point index for the bottom-side optical edge, or -1 if this glyph has none.",
            ),
        ],
    ),
    #
    # TSIC
    #
    (
        "TSIC",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version of table initially set to 0x00010000.",
            ),
            FieldSpec("uint16", "Flags", description="TSIC flags - set to 0"),
            FieldSpec("uint16", "AxisCount", description="Axis count from fvar"),
            FieldSpec("uint16", "RecordCount", description="TSIC record count"),
            FieldSpec("uint16", "Reserved", description="Set to 0"),
            FieldSpec(
                "Tag",
                "AxisArray",
                repeat="AxisCount",
                aux=0,
                description="Array of axis tags in fvar order",
            ),
            FieldSpec(
                "LocationRecord",
                "RecordLocations",
                repeat="RecordCount",
                aux=0,
                description="Location in variation space of TSIC record",
            ),
            FieldSpec(
                "TSICRecord",
                "Record",
                repeat="RecordCount",
                aux=0,
                description="Array of TSIC records",
            ),
        ],
    ),
    (
        "LocationRecord",
        [
            FieldSpec(
                "F2Dot14", "Axis", repeat="AxisCount", aux=0, description="Axis record"
            ),
        ],
    ),
    (
        "TSICRecord",
        [
            FieldSpec("uint16", "Flags", description="Record flags - set to 0"),
            FieldSpec(
                "uint16",
                "NumCVTEntries",
                description="Number of CVT number value pairs",
            ),
            FieldSpec(
                "uint16",
                "NameLength",
                description="Length of optional user record name",
            ),
            FieldSpec(
                "uint16",
                "NameArray",
                repeat="NameLength",
                aux=0,
                description="Unicode 16 name",
            ),
            FieldSpec(
                "uint16",
                "CVTArray",
                repeat="NumCVTEntries",
                aux=0,
                description="CVT number array",
            ),
            FieldSpec(
                "int16",
                "CVTValueArray",
                repeat="NumCVTEntries",
                aux=0,
                description="CVT value",
            ),
        ],
    ),
    #
    # COLR
    #
    (
        "COLR",
        [
            FieldSpec(
                "uint16", "Version", description="Table version number (starts at 0)."
            ),
            FieldSpec(
                "uint16",
                "BaseGlyphRecordCount",
                description="Number of Base Glyph Records.",
            ),
            FieldSpec(
                "LOffset",
                "BaseGlyphRecordArray",
                description="Offset (from beginning of COLR table) to Base Glyph records.",
            ),
            FieldSpec(
                "LOffset",
                "LayerRecordArray",
                description="Offset (from beginning of COLR table) to Layer Records.",
            ),
            FieldSpec(
                "uint16", "LayerRecordCount", description="Number of Layer Records."
            ),
            FieldSpec(
                "LOffset",
                "BaseGlyphList",
                aux="Version >= 1",
                description="Offset (from beginning of COLR table) to array of Version-1 Base Glyph records.",
            ),
            FieldSpec(
                "LOffset",
                "LayerList",
                aux="Version >= 1",
                description="Offset (from beginning of COLR table) to LayerList.",
            ),
            FieldSpec(
                "LOffset",
                "ClipList",
                aux="Version >= 1",
                description="Offset to ClipList table (may be NULL)",
            ),
            FieldSpec(
                "LOffsetTo(DeltaSetIndexMap)",
                "VarIndexMap",
                aux="Version >= 1",
                description="Offset to DeltaSetIndexMap table (may be NULL)",
            ),
            FieldSpec(
                "LOffset",
                "VarStore",
                aux="Version >= 1",
                description="Offset to variation store (may be NULL)",
            ),
        ],
    ),
    (
        "BaseGlyphRecordArray",
        [
            FieldSpec(
                "BaseGlyphRecord",
                "BaseGlyphRecord",
                repeat="BaseGlyphRecordCount",
                aux=0,
                description="Base Glyph records.",
            ),
        ],
    ),
    (
        "BaseGlyphRecord",
        [
            FieldSpec(
                "GlyphID",
                "BaseGlyph",
                description="Glyph ID of reference glyph. This glyph is for reference only and is not rendered for color.",
            ),
            FieldSpec(
                "uint16",
                "FirstLayerIndex",
                description="Index (from beginning of the Layer Records) to the layer record. There will be numLayers consecutive entries for this base glyph.",
            ),
            FieldSpec(
                "uint16",
                "NumLayers",
                description="Number of color layers associated with this glyph.",
            ),
        ],
    ),
    (
        "LayerRecordArray",
        [
            FieldSpec(
                "LayerRecord",
                "LayerRecord",
                repeat="LayerRecordCount",
                aux=0,
                description="Layer records.",
            ),
        ],
    ),
    (
        "LayerRecord",
        [
            FieldSpec(
                "GlyphID",
                "LayerGlyph",
                description="Glyph ID of layer glyph (must be in z-order from bottom to top).",
            ),
            FieldSpec(
                "uint16",
                "PaletteIndex",
                description="Index value to use with a selected color palette.",
            ),
        ],
    ),
    (
        "BaseGlyphList",
        [
            FieldSpec(
                "uint32",
                "BaseGlyphCount",
                description="Number of Version-1 Base Glyph records",
            ),
            FieldSpec(
                "struct",
                "BaseGlyphPaintRecord",
                repeat="BaseGlyphCount",
                aux=0,
                description="Array of Version-1 Base Glyph records",
            ),
        ],
    ),
    (
        "BaseGlyphPaintRecord",
        [
            FieldSpec(
                "GlyphID", "BaseGlyph", description="Glyph ID of reference glyph."
            ),
            FieldSpec(
                "LOffset",
                "Paint",
                description="Offset (from beginning of BaseGlyphPaintRecord) to Paint, typically a PaintColrLayers.",
            ),
        ],
    ),
    (
        "LayerList",
        [
            FieldSpec("uint32", "LayerCount", description="Number of Version-1 Layers"),
            FieldSpec(
                "LOffset",
                "Paint",
                repeat="LayerCount",
                aux=0,
                description="Array of offsets to Paint tables, from the start of the LayerList table.",
            ),
        ],
    ),
    (
        "ClipListFormat1",
        [
            FieldSpec(
                "uint8",
                "Format",
                description="Format for ClipList with 16bit glyph IDs: 1",
            ),
            FieldSpec("uint32", "ClipCount", description="Number of Clip records."),
            FieldSpec(
                "struct",
                "ClipRecord",
                repeat="ClipCount",
                aux=0,
                description="Array of Clip records sorted by glyph ID.",
            ),
        ],
    ),
    (
        "ClipRecord",
        [
            FieldSpec(
                "uint16", "StartGlyphID", description="First glyph ID in the range."
            ),
            FieldSpec(
                "uint16", "EndGlyphID", description="Last glyph ID in the range."
            ),
            FieldSpec("Offset24", "ClipBox", description="Offset to a ClipBox table."),
        ],
    ),
    (
        "ClipBoxFormat1",
        [
            FieldSpec(
                "uint8",
                "Format",
                description="Format for ClipBox without variation: set to 1.",
            ),
            FieldSpec("int16", "xMin", description="Minimum x of clip box."),
            FieldSpec("int16", "yMin", description="Minimum y of clip box."),
            FieldSpec("int16", "xMax", description="Maximum x of clip box."),
            FieldSpec("int16", "yMax", description="Maximum y of clip box."),
        ],
    ),
    (
        "ClipBoxFormat2",
        [
            FieldSpec(
                "uint8", "Format", description="Format for variable ClipBox: set to 2."
            ),
            FieldSpec(
                "int16", "xMin", description="Minimum x of clip box. VarIndexBase + 0."
            ),
            FieldSpec(
                "int16", "yMin", description="Minimum y of clip box. VarIndexBase + 1."
            ),
            FieldSpec(
                "int16", "xMax", description="Maximum x of clip box. VarIndexBase + 2."
            ),
            FieldSpec(
                "int16", "yMax", description="Maximum y of clip box. VarIndexBase + 3."
            ),
            FieldSpec(
                "VarIndex",
                "VarIndexBase",
                description="Base index into DeltaSetIndexMap.",
            ),
        ],
    ),
    # COLRv1 Affine2x3 uses the same column-major order to serialize a 2D
    # Affine Transformation as the one used by fontTools.misc.transform.
    # However, for historical reasons, the labels 'xy' and 'yx' are swapped.
    # Their fundamental meaning is the same though.
    # COLRv1 Affine2x3 follows the names found in FreeType and Cairo.
    # In all case, the second element in the 6-tuple correspond to the
    # y-part of the x basis vector, and the third to the x-part of the y
    # basis vector.
    # See https://github.com/googlefonts/colr-gradients-spec/pull/85
    (
        "Affine2x3",
        [
            FieldSpec("Fixed", "xx", description="x-part of x basis vector"),
            FieldSpec("Fixed", "yx", description="y-part of x basis vector"),
            FieldSpec("Fixed", "xy", description="x-part of y basis vector"),
            FieldSpec("Fixed", "yy", description="y-part of y basis vector"),
            FieldSpec("Fixed", "dx", description="Translation in x direction"),
            FieldSpec("Fixed", "dy", description="Translation in y direction"),
        ],
    ),
    (
        "VarAffine2x3",
        [
            FieldSpec(
                "Fixed", "xx", description="x-part of x basis vector. VarIndexBase + 0."
            ),
            FieldSpec(
                "Fixed", "yx", description="y-part of x basis vector. VarIndexBase + 1."
            ),
            FieldSpec(
                "Fixed", "xy", description="x-part of y basis vector. VarIndexBase + 2."
            ),
            FieldSpec(
                "Fixed", "yy", description="y-part of y basis vector. VarIndexBase + 3."
            ),
            FieldSpec(
                "Fixed",
                "dx",
                description="Translation in x direction. VarIndexBase + 4.",
            ),
            FieldSpec(
                "Fixed",
                "dy",
                description="Translation in y direction. VarIndexBase + 5.",
            ),
            FieldSpec(
                "VarIndex",
                "VarIndexBase",
                description="Base index into DeltaSetIndexMap.",
            ),
        ],
    ),
    (
        "ColorStop",
        [
            FieldSpec("F2Dot14", "StopOffset"),
            FieldSpec(
                "uint16", "PaletteIndex", description="Index for a CPAL palette entry."
            ),
            FieldSpec(
                "F2Dot14", "Alpha", description="Values outsided [0.,1.] reserved"
            ),
        ],
    ),
    (
        "VarColorStop",
        [
            FieldSpec("F2Dot14", "StopOffset", description="VarIndexBase + 0."),
            FieldSpec(
                "uint16", "PaletteIndex", description="Index for a CPAL palette entry."
            ),
            FieldSpec(
                "F2Dot14",
                "Alpha",
                description="Values outsided [0.,1.] reserved. VarIndexBase + 1.",
            ),
            FieldSpec(
                "VarIndex",
                "VarIndexBase",
                description="Base index into DeltaSetIndexMap.",
            ),
        ],
    ),
    (
        "ColorLine",
        [
            FieldSpec(
                "ExtendMode",
                "Extend",
                description="Enum {PAD = 0, REPEAT = 1, REFLECT = 2}",
            ),
            FieldSpec("uint16", "StopCount", description="Number of Color stops."),
            FieldSpec(
                "ColorStop",
                "ColorStop",
                repeat="StopCount",
                aux=0,
                description="Array of Color stops.",
            ),
        ],
    ),
    (
        "VarColorLine",
        [
            FieldSpec(
                "ExtendMode",
                "Extend",
                description="Enum {PAD = 0, REPEAT = 1, REFLECT = 2}",
            ),
            FieldSpec("uint16", "StopCount", description="Number of Color stops."),
            FieldSpec(
                "VarColorStop",
                "ColorStop",
                repeat="StopCount",
                aux=0,
                description="Array of Color stops.",
            ),
        ],
    ),
    # PaintColrLayers
    (
        "PaintFormat1",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 1"
            ),
            FieldSpec(
                "uint8",
                "NumLayers",
                description="Number of offsets to Paint to read from LayerList.",
            ),
            FieldSpec("uint32", "FirstLayerIndex", description="Index into LayerList."),
        ],
    ),
    # PaintSolid
    (
        "PaintFormat2",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 2"
            ),
            FieldSpec(
                "uint16", "PaletteIndex", description="Index for a CPAL palette entry."
            ),
            FieldSpec(
                "F2Dot14", "Alpha", description="Values outsided [0.,1.] reserved"
            ),
        ],
    ),
    # PaintVarSolid
    (
        "PaintFormat3",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 3"
            ),
            FieldSpec(
                "uint16", "PaletteIndex", description="Index for a CPAL palette entry."
            ),
            FieldSpec(
                "F2Dot14",
                "Alpha",
                description="Values outsided [0.,1.] reserved. VarIndexBase + 0.",
            ),
            FieldSpec(
                "VarIndex",
                "VarIndexBase",
                description="Base index into DeltaSetIndexMap.",
            ),
        ],
    ),
    # PaintLinearGradient
    (
        "PaintFormat4",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 4"
            ),
            FieldSpec(
                "Offset24",
                "ColorLine",
                description="Offset (from beginning of PaintLinearGradient table) to ColorLine subtable.",
            ),
            FieldSpec("int16", "x0"),
            FieldSpec("int16", "y0"),
            FieldSpec("int16", "x1"),
            FieldSpec("int16", "y1"),
            FieldSpec("int16", "x2"),
            FieldSpec("int16", "y2"),
        ],
    ),
    # PaintVarLinearGradient
    (
        "PaintFormat5",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 5"
            ),
            FieldSpec(
                "LOffset24To(VarColorLine)",
                "ColorLine",
                description="Offset (from beginning of PaintVarLinearGradient table) to VarColorLine subtable.",
            ),
            FieldSpec("int16", "x0", description="VarIndexBase + 0."),
            FieldSpec("int16", "y0", description="VarIndexBase + 1."),
            FieldSpec("int16", "x1", description="VarIndexBase + 2."),
            FieldSpec("int16", "y1", description="VarIndexBase + 3."),
            FieldSpec("int16", "x2", description="VarIndexBase + 4."),
            FieldSpec("int16", "y2", description="VarIndexBase + 5."),
            FieldSpec(
                "VarIndex",
                "VarIndexBase",
                description="Base index into DeltaSetIndexMap.",
            ),
        ],
    ),
    # PaintRadialGradient
    (
        "PaintFormat6",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 6"
            ),
            FieldSpec(
                "Offset24",
                "ColorLine",
                description="Offset (from beginning of PaintRadialGradient table) to ColorLine subtable.",
            ),
            FieldSpec("int16", "x0"),
            FieldSpec("int16", "y0"),
            FieldSpec("uint16", "r0"),
            FieldSpec("int16", "x1"),
            FieldSpec("int16", "y1"),
            FieldSpec("uint16", "r1"),
        ],
    ),
    # PaintVarRadialGradient
    (
        "PaintFormat7",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 7"
            ),
            FieldSpec(
                "LOffset24To(VarColorLine)",
                "ColorLine",
                description="Offset (from beginning of PaintVarRadialGradient table) to VarColorLine subtable.",
            ),
            FieldSpec("int16", "x0", description="VarIndexBase + 0."),
            FieldSpec("int16", "y0", description="VarIndexBase + 1."),
            FieldSpec("uint16", "r0", description="VarIndexBase + 2."),
            FieldSpec("int16", "x1", description="VarIndexBase + 3."),
            FieldSpec("int16", "y1", description="VarIndexBase + 4."),
            FieldSpec("uint16", "r1", description="VarIndexBase + 5."),
            FieldSpec(
                "VarIndex",
                "VarIndexBase",
                description="Base index into DeltaSetIndexMap.",
            ),
        ],
    ),
    # PaintSweepGradient
    (
        "PaintFormat8",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 8"
            ),
            FieldSpec(
                "Offset24",
                "ColorLine",
                description="Offset (from beginning of PaintSweepGradient table) to ColorLine subtable.",
            ),
            FieldSpec("int16", "centerX", description="Center x coordinate."),
            FieldSpec("int16", "centerY", description="Center y coordinate."),
            FieldSpec(
                "BiasedAngle",
                "startAngle",
                description="Start of the angular range of the gradient.",
            ),
            FieldSpec(
                "BiasedAngle",
                "endAngle",
                description="End of the angular range of the gradient.",
            ),
        ],
    ),
    # PaintVarSweepGradient
    (
        "PaintFormat9",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 9"
            ),
            FieldSpec(
                "LOffset24To(VarColorLine)",
                "ColorLine",
                description="Offset (from beginning of PaintVarSweepGradient table) to VarColorLine subtable.",
            ),
            FieldSpec(
                "int16", "centerX", description="Center x coordinate. VarIndexBase + 0."
            ),
            FieldSpec(
                "int16", "centerY", description="Center y coordinate. VarIndexBase + 1."
            ),
            FieldSpec(
                "BiasedAngle",
                "startAngle",
                description="Start of the angular range of the gradient. VarIndexBase + 2.",
            ),
            FieldSpec(
                "BiasedAngle",
                "endAngle",
                description="End of the angular range of the gradient. VarIndexBase + 3.",
            ),
            FieldSpec(
                "VarIndex",
                "VarIndexBase",
                description="Base index into DeltaSetIndexMap.",
            ),
        ],
    ),
    # PaintGlyph
    (
        "PaintFormat10",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 10"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintGlyph table) to Paint subtable.",
            ),
            FieldSpec(
                "GlyphID", "Glyph", description="Glyph ID for the source outline."
            ),
        ],
    ),
    # PaintColrGlyph
    (
        "PaintFormat11",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 11"
            ),
            FieldSpec(
                "GlyphID",
                "Glyph",
                description="Virtual glyph ID for a BaseGlyphList base glyph.",
            ),
        ],
    ),
    # PaintTransform
    (
        "PaintFormat12",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 12"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintTransform table) to Paint subtable.",
            ),
            FieldSpec(
                "LOffset24To(Affine2x3)",
                "Transform",
                description="2x3 matrix for 2D affine transformations.",
            ),
        ],
    ),
    # PaintVarTransform
    (
        "PaintFormat13",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 13"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintVarTransform table) to Paint subtable.",
            ),
            FieldSpec(
                "LOffset24To(VarAffine2x3)",
                "Transform",
                description="2x3 matrix for 2D affine transformations.",
            ),
        ],
    ),
    # PaintTranslate
    (
        "PaintFormat14",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 14"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintTranslate table) to Paint subtable.",
            ),
            FieldSpec("int16", "dx", description="Translation in x direction."),
            FieldSpec("int16", "dy", description="Translation in y direction."),
        ],
    ),
    # PaintVarTranslate
    (
        "PaintFormat15",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 15"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintVarTranslate table) to Paint subtable.",
            ),
            FieldSpec(
                "int16",
                "dx",
                description="Translation in x direction. VarIndexBase + 0.",
            ),
            FieldSpec(
                "int16",
                "dy",
                description="Translation in y direction. VarIndexBase + 1.",
            ),
            FieldSpec(
                "VarIndex",
                "VarIndexBase",
                description="Base index into DeltaSetIndexMap.",
            ),
        ],
    ),
    # PaintScale
    (
        "PaintFormat16",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 16"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintScale table) to Paint subtable.",
            ),
            FieldSpec("F2Dot14", "scaleX"),
            FieldSpec("F2Dot14", "scaleY"),
        ],
    ),
    # PaintVarScale
    (
        "PaintFormat17",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 17"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintVarScale table) to Paint subtable.",
            ),
            FieldSpec("F2Dot14", "scaleX", description="VarIndexBase + 0."),
            FieldSpec("F2Dot14", "scaleY", description="VarIndexBase + 1."),
            FieldSpec(
                "VarIndex",
                "VarIndexBase",
                description="Base index into DeltaSetIndexMap.",
            ),
        ],
    ),
    # PaintScaleAroundCenter
    (
        "PaintFormat18",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 18"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintScaleAroundCenter table) to Paint subtable.",
            ),
            FieldSpec("F2Dot14", "scaleX"),
            FieldSpec("F2Dot14", "scaleY"),
            FieldSpec("int16", "centerX"),
            FieldSpec("int16", "centerY"),
        ],
    ),
    # PaintVarScaleAroundCenter
    (
        "PaintFormat19",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 19"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintVarScaleAroundCenter table) to Paint subtable.",
            ),
            FieldSpec("F2Dot14", "scaleX", description="VarIndexBase + 0."),
            FieldSpec("F2Dot14", "scaleY", description="VarIndexBase + 1."),
            FieldSpec("int16", "centerX", description="VarIndexBase + 2."),
            FieldSpec("int16", "centerY", description="VarIndexBase + 3."),
            FieldSpec(
                "VarIndex",
                "VarIndexBase",
                description="Base index into DeltaSetIndexMap.",
            ),
        ],
    ),
    # PaintScaleUniform
    (
        "PaintFormat20",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 20"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintScaleUniform table) to Paint subtable.",
            ),
            FieldSpec("F2Dot14", "scale"),
        ],
    ),
    # PaintVarScaleUniform
    (
        "PaintFormat21",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 21"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintVarScaleUniform table) to Paint subtable.",
            ),
            FieldSpec("F2Dot14", "scale", description="VarIndexBase + 0."),
            FieldSpec(
                "VarIndex",
                "VarIndexBase",
                description="Base index into DeltaSetIndexMap.",
            ),
        ],
    ),
    # PaintScaleUniformAroundCenter
    (
        "PaintFormat22",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 22"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintScaleUniformAroundCenter table) to Paint subtable.",
            ),
            FieldSpec("F2Dot14", "scale"),
            FieldSpec("int16", "centerX"),
            FieldSpec("int16", "centerY"),
        ],
    ),
    # PaintVarScaleUniformAroundCenter
    (
        "PaintFormat23",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 23"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintVarScaleUniformAroundCenter table) to Paint subtable.",
            ),
            FieldSpec("F2Dot14", "scale", description="VarIndexBase + 0"),
            FieldSpec("int16", "centerX", description="VarIndexBase + 1"),
            FieldSpec("int16", "centerY", description="VarIndexBase + 2"),
            FieldSpec(
                "VarIndex",
                "VarIndexBase",
                description="Base index into DeltaSetIndexMap.",
            ),
        ],
    ),
    # PaintRotate
    (
        "PaintFormat24",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 24"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintRotate table) to Paint subtable.",
            ),
            FieldSpec("Angle", "angle"),
        ],
    ),
    # PaintVarRotate
    (
        "PaintFormat25",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 25"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintVarRotate table) to Paint subtable.",
            ),
            FieldSpec("Angle", "angle", description="VarIndexBase + 0."),
            FieldSpec(
                "VarIndex",
                "VarIndexBase",
                description="Base index into DeltaSetIndexMap.",
            ),
        ],
    ),
    # PaintRotateAroundCenter
    (
        "PaintFormat26",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 26"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintRotateAroundCenter table) to Paint subtable.",
            ),
            FieldSpec("Angle", "angle"),
            FieldSpec("int16", "centerX"),
            FieldSpec("int16", "centerY"),
        ],
    ),
    # PaintVarRotateAroundCenter
    (
        "PaintFormat27",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 27"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintVarRotateAroundCenter table) to Paint subtable.",
            ),
            FieldSpec("Angle", "angle", description="VarIndexBase + 0."),
            FieldSpec("int16", "centerX", description="VarIndexBase + 1."),
            FieldSpec("int16", "centerY", description="VarIndexBase + 2."),
            FieldSpec(
                "VarIndex",
                "VarIndexBase",
                description="Base index into DeltaSetIndexMap.",
            ),
        ],
    ),
    # PaintSkew
    (
        "PaintFormat28",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 28"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintSkew table) to Paint subtable.",
            ),
            FieldSpec("Angle", "xSkewAngle"),
            FieldSpec("Angle", "ySkewAngle"),
        ],
    ),
    # PaintVarSkew
    (
        "PaintFormat29",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 29"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintVarSkew table) to Paint subtable.",
            ),
            FieldSpec("Angle", "xSkewAngle", description="VarIndexBase + 0."),
            FieldSpec("Angle", "ySkewAngle", description="VarIndexBase + 1."),
            FieldSpec(
                "VarIndex",
                "VarIndexBase",
                description="Base index into DeltaSetIndexMap.",
            ),
        ],
    ),
    # PaintSkewAroundCenter
    (
        "PaintFormat30",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 30"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintSkewAroundCenter table) to Paint subtable.",
            ),
            FieldSpec("Angle", "xSkewAngle"),
            FieldSpec("Angle", "ySkewAngle"),
            FieldSpec("int16", "centerX"),
            FieldSpec("int16", "centerY"),
        ],
    ),
    # PaintVarSkewAroundCenter
    (
        "PaintFormat31",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 31"
            ),
            FieldSpec(
                "Offset24",
                "Paint",
                description="Offset (from beginning of PaintVarSkewAroundCenter table) to Paint subtable.",
            ),
            FieldSpec("Angle", "xSkewAngle", description="VarIndexBase + 0."),
            FieldSpec("Angle", "ySkewAngle", description="VarIndexBase + 1."),
            FieldSpec("int16", "centerX", description="VarIndexBase + 2."),
            FieldSpec("int16", "centerY", description="VarIndexBase + 3."),
            FieldSpec(
                "VarIndex",
                "VarIndexBase",
                description="Base index into DeltaSetIndexMap.",
            ),
        ],
    ),
    # PaintComposite
    (
        "PaintFormat32",
        [
            FieldSpec(
                "uint8", "PaintFormat", description="Format identifier-format = 32"
            ),
            FieldSpec(
                "LOffset24To(Paint)",
                "SourcePaint",
                description="Offset (from beginning of PaintComposite table) to source Paint subtable.",
            ),
            FieldSpec(
                "CompositeMode",
                "CompositeMode",
                description="A CompositeMode enumeration value.",
            ),
            FieldSpec(
                "LOffset24To(Paint)",
                "BackdropPaint",
                description="Offset (from beginning of PaintComposite table) to backdrop Paint subtable.",
            ),
        ],
    ),
    #
    # avar
    #
    (
        "AxisValueMap",
        [
            FieldSpec(
                "F2Dot14",
                "FromCoordinate",
                description="A normalized coordinate value obtained using default normalization",
            ),
            FieldSpec(
                "F2Dot14",
                "ToCoordinate",
                description="The modified, normalized coordinate value",
            ),
        ],
    ),
    (
        "AxisSegmentMap",
        [
            FieldSpec(
                "uint16",
                "PositionMapCount",
                description="The number of correspondence pairs for this axis",
            ),
            FieldSpec(
                "AxisValueMap",
                "AxisValueMap",
                repeat="PositionMapCount",
                aux=0,
                description="The array of axis value map records for this axis",
            ),
        ],
    ),
    (
        "avar",
        [
            FieldSpec(
                "Version",
                "Version",
                description="Version of the avar table- 0x00010000 or 0x00020000",
            ),
            FieldSpec(
                "uint16", "Reserved", description="Permanently reserved; set to zero"
            ),
            FieldSpec(
                "uint16",
                "AxisCount",
                description='The number of variation axes for this font. This must be the same number as axisCount in the "fvar" table',
            ),
            FieldSpec(
                "AxisSegmentMap",
                "AxisSegmentMap",
                repeat="AxisCount",
                aux=0,
                description='The segment maps array — one segment map for each axis, in the order of axes specified in the "fvar" table',
            ),
            FieldSpec(
                "LOffsetTo(DeltaSetIndexMap)", "VarIdxMap", aux="Version >= 0x00020000"
            ),
            FieldSpec("LOffset", "VarStore", aux="Version >= 0x00020000"),
        ],
    ),
    #
    # IFT - Incremental Font Transfer tables
    # https://w3c.github.io/IFT/Overview.html
    # Reference: https://github.com/googlefonts/fontations/blob/main/read-fonts/src/tables/ift.rs
    #
    (
        "PatchMapFormat2",
        [
            FieldSpec(
                "uint8", "Format", description="Set to 2, identifies this as format 2."
            ),
            FieldSpec("uint24", "Reserved", description="Not used, set to 0."),
            FieldSpec(
                "uint8",
                "Flags",
                description="Bit 0: CffCharStringsOffset present. Bit 1: Cff2CharStringsOffset present.",
            ),
            FieldSpec(
                "uint32",
                "CompatibilityId",
                repeat=4,
                aux=0,
                description="Unique ID to identify compatible patches (16 bytes).",
            ),
            FieldSpec(
                "uint8",
                "DefaultPatchFormat",
                description="Default format of patches linked to by urlTemplate.",
            ),
            FieldSpec(
                "uint24",
                "NumEntries",
                description="Number of entries encoded in this table.",
            ),
            FieldSpec(
                "LOffset",
                "MappingEntries",
                description="Offset to a MappingEntries sub-table.",
            ),
            FieldSpec(
                "LOffset",
                "EntryIdStringData",
                description="Offset to entry ID string data block. May be null (0).",
            ),
            FieldSpec(
                "uint16",
                "UrlTemplateLength",
                description="Length of the urlTemplate byte array.",
            ),
            FieldSpec(
                "uint8",
                "UrlTemplate",
                repeat="UrlTemplateLength",
                aux=0,
                description="URL Template bytes used to produce URL strings for each entry.",
            ),
            FieldSpec(
                "uint32",
                "CffCharStringsOffset",
                aux="Flags & 0x01",
                description="Offset from start of CFF table to CharStrings INDEX.",
            ),
            FieldSpec(
                "uint32",
                "Cff2CharStringsOffset",
                aux="Flags & 0x02",
                description="Offset from start of CFF2 table to CharStrings INDEX.",
            ),
        ],
    ),
    # 'MappingEntries' contains stateful, delta-encoded entry IDs.  The custom
    # 'MappingEntriesConverter' is used to do a stateful parsing of all records.
    (
        "MappingEntries",
        [
            FieldSpec(
                "MappingEntries",
                "entries",
                description="Array of MappingEntry records.",
            )
        ],
    ),
    ("EntryIdStringData", []),
]
