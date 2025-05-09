# Merge Log.final.out files from multiple runs into one table
# usage:
# awk -f merge_star_out.awk */Log.final.out

BEGIN {
    FS="|";
    for (jj=1;jj<=ARGC;jj++)
    {
        a=ARGV[jj];
        gsub("/Log.final.out","",a);
        printf ";" a
    };
    printf "\n";
}
{
    gsub(/^[ \t]+|[ \t]+$/,"",$1);
    gsub(/^[ \t]+|[ \t]+$/,"",$2);
    L[FNR]=$1;
    V[FNR,ARGIND]=$2
}
END {
    for (ii=1;ii<=length(L);ii++)
    {
        printf "%s",L[ii];
        if (V[ii,1]!="")
            for (jj=1;jj<=ARGC;jj++)
                printf ";" V[ii,jj];
        printf "\n"
     }
}
