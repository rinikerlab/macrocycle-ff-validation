#!/bin/bash

AWK_PARSE_FLAG='
BEGIN{
    current_flag="none"
}
NF==3 && $1 == "[" && $3 == "]" {
    current_flag=$2
}
'

AWK_REMOVE="
$AWK_PARSE_FLAG
current_flag != flag{print}
"

AWK_PREPEND="
$AWK_PARSE_FLAG
(current_flag == flag) && (previous != flag) {
    printf \"%s\\n\",string;
}
{
    previous=current_flag;
    print;
}
"

AWK_APPEND="
$AWK_PARSE_FLAG
(current_flag != flag) && (previous == flag) {
    printf \"%s\\n\",string;
}
{
    previous=current_flag;
    print;
}
"

AWK_EXTRACT="
$AWK_PARSE_FLAG
(current_flag == flag) && (\$1 != \"[\") {print}
"

if [ "$1" == "remove" ]; then
    if [ -z "$2" ]; then
        echo "usage: $0 remove FLAG [TOP]"
        exit 1
    fi
    flag="$2"
    awk -v flag=$flag "$AWK_REMOVE" $3
elif [ "$1" == "extract" ]; then
    if [ -z "$2" ]; then
        echo "usage: $0 extract FLAG [TOP]"
        exit 1
    fi
    flag="$2"
    awk -v flag=$flag "$AWK_EXTRACT" $3
elif [ "$1" == "prepend" ]; then
    if [ -z "$3" ]; then
        echo "usage: $0 prepend FLAG FSTRING [TOP]"
        exit 1
    fi
    flag="$2"
    string="$3"
    awk -v flag=$flag -v string="$string" "$AWK_PREPEND" $4
elif [ "$1" == "append" ]; then
    if [ -z "$3" ]; then
        echo "usage: $0 append FLAG FSTRING [TOP]"
        exit 1
    fi
    flag="$2"
    awk -v flag=$flag -v string="$3" "$AWK_APPEND" $4
else
    echo "commands: remove prepend; For help: $0 COMMAND"
fi
