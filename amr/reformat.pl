$grp = 0;
$mod = 8;

while(<>) {
    for(split(",")) {
        /(-?\d*)\.(\d*)F/ or next;
        push @ints, $1;
        push @fracs, $2;
    }
}

$ilen = (sort {$b <=> $a} map {length} @ints )[0];
$flen = (sort {$b <=> $a} map {length} @fracs)[0];

for($t = 0; $t <= $#fracs; $t++){
    if($grp>0 && !($t%$grp)) {
        printf "{";
    }
    printf "%${ilen}s.%-${flen}s", $ints[$t], $fracs[$t];
    if($grp>0 && $t!=0 && !(($t+1)%$grp)) {
        printf "}";
    }
    if(!(($t+1)%$mod)) {
        printf ",\n";
    }else {
        printf ", ";
    }
}
printf "\n";
