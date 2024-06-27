void frag() {
    vec4 mainCol = texture2D_0(vec2(TexCoord.x, TexCoord.y));
    gl_FragColor = vec4(1.0 - mainCol.r, 1.0 - mainCol.g, 1.0 - mainCol.b, mainCol.a);
}
