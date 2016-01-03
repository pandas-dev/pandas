// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_TYPES_H
#define PANDAS_TYPES_H

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace pandas {

enum TypeEnum {
  // A degerate NULL type
  NA = 0,

  // Little-endian integer types
  UINT8 = 1,
  INT8 = 2,
  UINT16 = 3,
  INT16 = 4,
  UINT32 = 5,
  INT32 = 6,
  UINT64 = 7,
  INT64 = 8,

  // A boolean value represented as 1 byte
  BOOL = 9,

  // 4-byte floating point value
  FLOAT = 10,

  // 8-byte floating point value
  DOUBLE = 11,

  // PyObject*
  PYOBJECT = 12,

  // Timestamp types
  TIMESTAMP = 13,
  TIMESTAMP_TZ = 14,

  // UTF8 variable-length string
  STRING = 15,

  // Categorical
  CATEGORY = 16
};


struct DataType {

  TypeEnum type;

  explicit DataType(TypeEnum type)
      : type(type) {}

  virtual std::string ToString() const = 0;

  virtual bool Equals(const DataType& other) {
    return type == other.type;
  }
};


typedef std::shared_ptr<DataType> TypePtr;


struct TimestampType : public DataType {

  enum class Unit: char {
    NANO = 0
  };

  Unit unit;

  explicit TimestampType(Unit unit = Unit::NANO)
      : DataType(TypeEnum::TIMESTAMP),
        unit(unit) {}

  TimestampType(const TimestampType& other)
      : TimestampType(other.unit) {}

  static char const *name() {
    return "timestamp";
  }

  virtual std::string ToString() const {
    return name();
  }
};


struct PyObjectType : public DataType {

  explicit PyObjectType()
      : DataType(TypeEnum::PYOBJECT) {}

  PyObjectType(const PyObjectType& other)
      : PyObjectType() {}

  static char const *name() {
    return "object";
  }

  virtual std::string ToString() const {
     return name();
  }
};


template <typename Derived>
struct PrimitiveType : public DataType {
  PrimitiveType()
      : DataType(Derived::type_enum) {}

  virtual std::string ToString() const {
    return std::string(static_cast<const Derived*>(this)->name());
  }
};


#define PRIMITIVE_DECL(TYPENAME, C_TYPE, ENUM, SIZE, NAME)  \
  typedef C_TYPE c_type;                                    \
  static constexpr TypeEnum type_enum = TypeEnum::ENUM;     \
  static constexpr size_t size = SIZE;                      \
                                                            \
  explicit TYPENAME()                                       \
      : PrimitiveType<TYPENAME>() {}                        \
                                                            \
  static const char* name() {                               \
    return NAME;                                            \
  }


struct NullType : public PrimitiveType<NullType> {
  PRIMITIVE_DECL(NullType, void, NA, 0, "null");
};

struct UInt8Type : public PrimitiveType<UInt8Type> {
  PRIMITIVE_DECL(UInt8Type, uint8_t, UINT8, 1, "uint8");
};

struct Int8Type : public PrimitiveType<Int8Type> {
  PRIMITIVE_DECL(Int8Type, int8_t, INT8, 1, "int8");
};

struct UInt16Type : public PrimitiveType<UInt16Type> {
  PRIMITIVE_DECL(UInt16Type, uint16_t, UINT16, 2, "uint16");
};

struct Int16Type : public PrimitiveType<Int16Type> {
  PRIMITIVE_DECL(Int16Type, int16_t, INT16, 2, "int16");
};

struct UInt32Type : public PrimitiveType<UInt32Type> {
  PRIMITIVE_DECL(UInt32Type, uint32_t, UINT32, 4, "uint32");
};

struct Int32Type : public PrimitiveType<Int32Type> {
  PRIMITIVE_DECL(Int32Type, int32_t, INT32, 4, "int32");
};

struct UInt64Type : public PrimitiveType<UInt64Type> {
  PRIMITIVE_DECL(UInt64Type, uint64_t, UINT64, 8, "uint64");
};

struct Int64Type : public PrimitiveType<Int64Type> {
  PRIMITIVE_DECL(Int64Type, int64_t, INT64, 8, "int64");
};

struct FloatType : public PrimitiveType<FloatType> {
  PRIMITIVE_DECL(FloatType, float, FLOAT, 4, "float");
};

struct DoubleType : public PrimitiveType<DoubleType> {
  PRIMITIVE_DECL(DoubleType, double, DOUBLE, 8, "double");
};

struct BooleanType : public PrimitiveType<BooleanType> {
  PRIMITIVE_DECL(BooleanType, uint8_t, BOOL, 1, "bool");
};

} // namespace pandas

#endif  // PANDAS_TYPES_H
