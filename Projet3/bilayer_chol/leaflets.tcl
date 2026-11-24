proc load_leaflets {fname {molid top}} {
    if {$molid eq "top"} { set molid [molinfo top] }
    set fp [open $fname r]; set data [read $fp]; close $fp
    set cur ""
    foreach line [split $data "\n"] {
        set line [string trim $line]
        if {$line eq ""} continue
        if {[regexp {^\[\s*(.*?)\s*\]$} $line -> name]} {
            set cur $name; lappend order $name; set groups($cur) {}; continue
        }
        foreach x $line { lappend groups($cur) [expr {$x - 1}] }   ;# 1-based -> 0-based
    }
    set box [lindex [pbc get -molid $molid] 0]
    set newids {}
    foreach name $order {
        set sel [atomselect $molid "residue $groups($name)"]
        set tmp "tmp_$name.pdb"
        $sel writepdb $tmp
        puts "$name : [llength $groups($name)] lipids, [$sel num] atoms"
        $sel delete
        set new [mol new $tmp type pdb waitfor all]
        mol rename $new $name
        pbc set $box -molid $new
        file delete $tmp
        lappend newids $new
    }
    return $newids
}
