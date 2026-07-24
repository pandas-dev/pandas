/****************************************************************************
**
** Copyright (C) 2019 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

import QtQuick 2.15
import QtQuick3D 1.15
import QtQuick3D.Materials 1.15

CustomMaterial {
    // These properties names need to match the ones in the shader code!
    property real roughness: 0.0
    property real blur_size: 8.0
    property real refract_depth: 5
    property bool uEnvironmentMappingEnabled: true
    property bool uShadowMappingEnabled: false
    property real glass_bfactor: 0.0
    property bool glass_binside: false
    property real uFresnelPower: 1.0
    property real reflectivity_amount: 1.0
    property real glass_ior: 1.5
    property real intLightFall: 2.0
    property real intLightRot: 0.0
    property real intLightBrt: 0.0
    property real bumpScale: 0.5
    property int bumpBands: 1
    property vector3d bumpCoords: Qt.vector3d(1.0, 1.0, 1.0)
    property vector2d intLightPos: Qt.vector2d(0.5, 0.0)
    property vector3d glass_color: Qt.vector3d(0.9, 0.9, 0.9)
    property vector3d intLightCol: Qt.vector3d(0.9, 0.9, 0.9)
    hasTransparency: true

    shaderInfo: ShaderInfo {
        version: "330"
        type: "GLSL"
        shaderKey: ShaderInfo.Refraction | ShaderInfo.Glossy
    }

    property TextureInput glass_bump: TextureInput {
        enabled: true
        texture: Texture {
            id: glassBumpMap
            source: "maps/spherical_checker.png"
        }
    }

    property TextureInput uEnvironmentTexture: TextureInput {
            enabled: uEnvironmentMappingEnabled
            texture: Texture {
                id: envImage
                source: "maps/spherical_checker.png"
            }
    }
    property TextureInput uBakedShadowTexture: TextureInput {
            enabled: uShadowMappingEnabled
            texture: Texture {
                id: shadowImage
                source: "maps/shadow.png"
            }
    }
    property TextureInput randomGradient1D: TextureInput {
            texture: Texture {
                tilingModeHorizontal: Texture.Repeat
                tilingModeVertical: Texture.Repeat
                source: "maps/randomGradient1D.png"
            }
    }
    property TextureInput randomGradient2D: TextureInput {
            texture: Texture {
                tilingModeHorizontal: Texture.Repeat
                tilingModeVertical: Texture.Repeat
                source: "maps/randomGradient2D.png"
            }
    }
    property TextureInput randomGradient3D: TextureInput {
        texture: Texture {
            tilingModeHorizontal: Texture.Repeat
            tilingModeVertical: Texture.Repeat
            source: "maps/randomGradient3D.png"
        }
    }
    property TextureInput randomGradient4D: TextureInput {
        texture: Texture {
            tilingModeHorizontal: Texture.Repeat
            tilingModeVertical: Texture.Repeat
            source: "maps/randomGradient4D.png"
        }
    }

    Shader {
        id: mainShader
        stage: Shader.Fragment
        shader: "shaders/frostedThinGlass.frag"
    }
    Shader {
        id: noopShader
        stage: Shader.Fragment
        shader: "shaders/frostedThinGlassNoop.frag"
    }
    Shader {
        id: preBlurShader
        stage: Shader.Fragment
        shader: "shaders/frostedThinGlassPreBlur.frag"
    }
    Shader {
        id: blurXShader
        stage: Shader.Fragment
        shader: "shaders/frostedThinGlassBlurX.frag"
    }
    Shader {
        id: blurYShader
        stage: Shader.Fragment
        shader: "shaders/frostedThinGlassBlurY.frag"
    }

    Buffer {
        id: frameBuffer
        name: "frameBuffer"
        format: Buffer.Unknown
        textureFilterOperation: Buffer.Linear
        textureCoordOperation: Buffer.ClampToEdge
        sizeMultiplier: 1.0
        bufferFlags: Buffer.None // aka frame
    }

    Buffer {
        id: dummyBuffer
        name: "dummyBuffer"
        format: Buffer.RGBA8
        textureFilterOperation: Buffer.Linear
        textureCoordOperation: Buffer.ClampToEdge
        sizeMultiplier: 1.0
        bufferFlags: Buffer.None // aka frame
    }

    Buffer {
        id: tempBuffer
        name: "tempBuffer"
        format: Buffer.RGBA16F
        textureFilterOperation: Buffer.Linear
        textureCoordOperation: Buffer.ClampToEdge
        sizeMultiplier: 0.5
        bufferFlags: Buffer.None // aka frame
    }

    Buffer {
        id: blurYBuffer
        name: "tempBlurY"
        format: Buffer.RGBA16F
        textureFilterOperation: Buffer.Linear
        textureCoordOperation: Buffer.ClampToEdge
        sizeMultiplier: 0.5
        bufferFlags: Buffer.None // aka frame
    }

    Buffer {
        id: blurXBuffer
        name: "tempBlurX"
        format: Buffer.RGBA16F
        textureFilterOperation: Buffer.Linear
        textureCoordOperation: Buffer.ClampToEdge
        sizeMultiplier: 0.5
        bufferFlags: Buffer.None // aka frame
    }

    passes: [ Pass {
            shaders: noopShader
            output: dummyBuffer
            commands: [ BufferBlit {
                    destination: frameBuffer
                }
            ]
        }, Pass {
            shaders: preBlurShader
            output: tempBuffer
            commands: [ BufferInput {
                    buffer: frameBuffer
                    param: "OriginBuffer"
                }
            ]
        }, Pass {
            shaders: blurXShader
            output: blurXBuffer
            commands: [ BufferInput {
                    buffer: tempBuffer
                    param: "BlurBuffer"
                }
            ]
        }, Pass {
            shaders: blurYShader
            output: blurYBuffer
            commands: [ BufferInput {
                    buffer: blurXBuffer
                    param: "BlurBuffer"
                }, BufferInput {
                    buffer: tempBuffer
                    param: "OriginBuffer"
                }
            ]
        }, Pass {
            shaders: mainShader
            commands: [BufferInput {
                    buffer: blurYBuffer
                    param: "refractiveTexture"
                }, Blending {
                    srcBlending: Blending.SrcAlpha
                    destBlending: Blending.OneMinusSrcAlpha
                }
            ]
        }
    ]
}
