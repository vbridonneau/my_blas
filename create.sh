if [ $# -lt 1 ]
then
    echo "Usage $0 <routine_name>"
    exit 0
elif [ -f $1.h ]
then
    echo "$1.h already exist"
    exit 0
fi


routine=$1
upper=`echo $routine | tr [:lower:] [:upper:]`
# HEADER
echo -e "#ifndef DEF_${upper}_H" >  $routine.h
echo -e "#define DEF_${upper}_H" >> $routine.h
echo -e "\nvoid $routine();\n" >> $routine.h
echo -e "#endif//DEF_${upper}_H" >> $routine.h
# CODE
echo -e "#include \"$routine.h\"" > $routine.c
echo -e "\nvoid $routine() {\n\n}\n" >> $routine.c
# TEST FILE
echo -e "#include \"$routine.h\"" > $routine-perf.c
echo -e "\nvoid test_result() {\n\n}\n" >> $routine-perf.c
echo -e "int main(int argc, char **argv) {\n\ttest_result();" >> $routine-perf.c
echo -e "\treturn 0;\n}" >> $routine-perf.c
