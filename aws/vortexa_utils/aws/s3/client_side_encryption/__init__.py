# @Author: richard
# @Date:   2018-11-28T15:15:44+00:00
# @Last modified by:   richard
# @Last modified time: 2018-11-28T15:15:44+00:00

"""
# From the RUBY Docs.

Provides an encryption client that encrypts and decrypts data client-side,
storing the encrypted data in Amazon S3.

This client uses a process called "envelope encryption". Your private
encryption keys and your data's plain-text are **never** sent to
Amazon S3. **If you lose you encryption keys, you will not be able to
decrypt your data.**

## Envelope Encryption Overview

The goal of envelope encryption is to combine the performance of
fast symmetric encryption while maintaining the secure key management
that asymmetric keys provide.

A one-time-use symmetric key (envelope key) is generated client-side.
This is used to encrypt the data client-side. This key is then
encrypted by your master key and stored alongside your data in Amazon
S3.

When accessing your encrypted data with the encryption client,
the encrypted envelope key is retrieved and decrypted client-side
with your master key. The envelope key is then used to decrypt the
data client-side.

One of the benefits of envelope encryption is that if your master key
is compromised, you have the option of just re-encrypting the stored
envelope symmetric keys, instead of re-encrypting all of the
data in your account.

## Basic Usage

The encryption client requires an {Aws::S3::Client}. If you do not
provide a `:client`, then a client will be constructed for you.

    require 'openssl'
    key = OpenSSL::PKey::RSA.new(1024)

    # encryption client
    s3 = aws.s3.client_side_encryption.Client(encryption_key: key)

    # round-trip an object, encrypted/decrypted locally
    s3.put_object(bucket:'aws-sdk', key:'secret', body:'handshake')
    s3.get_object(bucket:'aws-sdk', key:'secret').body.read
    #=> 'handshake'

    # reading encrypted object without the encryption client
    # results in the getting the cipher text
    Aws::S3::Client.new.get_object(bucket:'aws-sdk', key:'secret').body.read
    #=> "... cipher text ..."

## Keys

For client-side encryption to work, you must provide one of the following:

* An encryption key
* A {KeyProvider}
* A KMS encryption key id

### An Encryption Key

You can pass a single encryption key. This is used as a master key
encrypting and decrypting all object keys.

    key = OpenSSL::Cipher.new("AES-256-ECB").random_key # symmetric key
    key = OpenSSL::PKey::RSA.new(1024) # asymmetric key pair

    s3 = Aws::S3::Encryption::Client.new(encryption_key: key)

### Key Provider

Alternatively, you can use a {KeyProvider}. A key provider makes
it easy to work with multiple keys and simplifies key rotation.

### KMS Encryption Key Id

If you pass the id to an AWS Key Management Service (KMS) key,
then KMS will be used to generate, encrypt and decrypt object keys.

    # keep track of the kms key id
    kms = Aws::KMS::Client.new
    key_id = kms.create_key.key_metadata.key_id

    Aws::S3::Encryption::Client.new(
      kms_key_id: key_id,
      kms_client: kms,
    )

## Custom Key Providers

A {KeyProvider} is any object that responds to:

* `#encryption_materials`
* `#key_for(materials_description)`

Here is a trivial implementation of an in-memory key provider.
This is provided as a demonstration of the key provider interface,
and should not be used in production:

    class KeyProvider

      def initialize(default_key_name, keys)
        @keys = keys
        @encryption_materials = Aws::S3::Encryption::Materials.new(
          key: @keys[default_key_name],
          description: JSON.dump(key: default_key_name),
        )
      end

      attr_reader :encryption_materials

      def key_for(matdesc)
        key_name = JSON.load(matdesc)['key']
        if key = @keys[key_name]
          key
        else
          raise "encryption key not found for: #{matdesc.inspect}"
        end
      end
    end

Given the above key provider, you can create an encryption client that
chooses the key to use based on the materials description stored with
the encrypted object. This makes it possible to use multiple keys
and simplifies key rotation.

    # uses "new-key" for encrypting objects, uses either for decrypting
    keys = KeyProvider.new('new-key', {
      "old-key" => Base64.decode64("kM5UVbhE/4rtMZJfsadYEdm2vaKFsmV2f5+URSeUCV4="),
      "new-key" => Base64.decode64("w1WLio3agRWRTSJK/Ouh8NHoqRQ6fn5WbSXDTHjXMSo="),
    }),

    # chooses the key based on the materials description stored
    # with the encrypted object
    s3 = Aws::S3::Encryption::Client.new(key_provider: keys)

## Materials Description

A materials description is JSON document string that is stored
in the metadata (or instruction file) of an encrypted object.
The {DefaultKeyProvider} uses the empty JSON document `"{}"`.

When building a key provider, you are free to store whatever
information you need to identify the master key that was used
to encrypt the object.

## Envelope Location

By default, the encryption client store the encryption envelope
with the object, as metadata. You can choose to have the envelope
stored in a separate "instruction file". An instruction file
is an object, with the key of the encrypted object, suffixed with
`".instruction"`.

Specify the `:envelope_location` option as `:instruction_file` to
use an instruction file for storing the envelope.

    # default behavior
    s3 = Aws::S3::Encryption::Client.new(
      key_provider: ...,
      envelope_location: :metadata,
    )

    # store envelope in a separate object
    s3 = Aws::S3::Encryption::Client.new(
      key_provider: ...,
      envelope_location: :instruction_file,
      instruction_file_suffix: '.instruction' # default
    )

When using an instruction file, multiple requests are made when
putting and getting the object. **This may cause issues if you are
issuing concurrent PUT and GET requests to an encrypted object.**
"""

from .client import Client
