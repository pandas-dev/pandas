-- see dviread._LuatexKpsewhich
kpse.set_program_name("latex")
kpse.init_prog("", 600, "ljfour")
while true do print(kpse.lookup(io.read():gsub("\r", ""), {mktexpk=true})); io.flush(); end
