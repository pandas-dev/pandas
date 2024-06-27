/****************************************************************************
**
** Copyright (C) 2020 The Qt Company Ltd.
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
import QtQuick3D.Effects 1.15

Effect {
    readonly property TextureInput downsample2: TextureInput {
        texture: Texture {}
    }
    readonly property TextureInput downsample4: TextureInput {
        texture: Texture {}
    }
    property real gamma: 1              // 0.1 - 4
    property real exposure: 0           // -9 - 9
    readonly property real exposureExp2: Math.pow(2, exposure)
    property real bloomThreshold: 1
    property real blurFalloff: 0        // 0 - 10
    readonly property real negativeBlurFalloffExp2: Math.pow(2, -blurFalloff)
    property real tonemappingLerp: 1    // 0 - 1
    property real channelThreshold: 1
    readonly property real poissonRotation: 0
    readonly property real poissonDistance: 4

    Shader {
        id: luminosityVert
        stage: Shader.Vertex
        shader: "shaders/luminosity.vert"
    }
    Shader {
        id: luminosityFrag
        stage: Shader.Fragment
        shader: "shaders/luminosity.frag"
    }

    Shader {
        id: blurVert
        stage: Shader.Vertex
        shader: "shaders/poissonblur.vert"
    }
    Shader {
        id: blurFrag
        stage: Shader.Fragment
        shader: "shaders/poissonblur.frag"
    }

    Shader {
        id: combiner
        stage: Shader.Fragment
        shader: "shaders/combiner.frag"
    }

    Buffer {
        id: luminosity_buffer2
        name: "luminosity_buffer2"
        format: Buffer.RGBA8
        textureFilterOperation: Buffer.Linear
        textureCoordOperation: Buffer.ClampToEdge
        bufferFlags: Buffer.None
        sizeMultiplier: 0.5
    }
    Buffer {
        id: downsample_buffer2
        name: "downsample_buffer2"
        format: Buffer.RGBA8
        textureFilterOperation: Buffer.Linear
        textureCoordOperation: Buffer.ClampToEdge
        bufferFlags: Buffer.None
        sizeMultiplier: 0.5
    }
    Buffer {
        id: downsample_buffer4
        name: "downsample_buffer4"
        format: Buffer.RGBA8
        textureFilterOperation: Buffer.Linear
        textureCoordOperation: Buffer.ClampToEdge
        bufferFlags: Buffer.None
        sizeMultiplier: 0.25
    }

    passes: [
        Pass {
            shaders: [ luminosityVert, luminosityFrag ]
            output: downsample_buffer2
        },
        Pass {
            shaders: [ luminosityVert, luminosityFrag ]
            commands: BufferInput {
                buffer: downsample_buffer2
            }
            output: luminosity_buffer2
        },
        Pass {
            shaders: [ blurVert, blurFrag ]
            commands: BufferInput {
                buffer: luminosity_buffer2
            }
            output: downsample_buffer2
        },
        Pass {

            shaders: [ blurVert, blurFrag ]
            commands: [
                SetUniformValue {
                    target: "poissonRotation"
                    value: 0.62831
                },
                BufferInput {
                    buffer: luminosity_buffer2
                }
            ]
            output: downsample_buffer4
        },
        Pass {
            shaders: combiner
            commands: [
                BufferInput {
                    param: "downsample2"
                    buffer: downsample_buffer2
                },
                BufferInput {
                    param: "downsample4"
                    buffer: downsample_buffer4
                }
            ]
        }
    ]
}
