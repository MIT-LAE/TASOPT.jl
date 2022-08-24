" simple comands to speed up fortran 2 julia conversion

"replace comment characters
:%s/^[c]/#/gi
:%s/!/#/g

"Replace subroutines with functions
:%s/subroutine/function/gi

"Replace & line continuations with just a continuation
:%s/\n\s*&\s*/\r\t/g

"Move any operator to previous line
:%s/\n\s*&\s*+/ +\r\t/g
:%s/\n\s*+/ +\r\t/g

:%s/\n\s*&\s*\*/ \*\r\t/g
:%s/\n\s*\*/ \*\r\t/g

:%s/\n\s*&\s*\// \/\r\t/g
:%s/\n\s*\// \/\r\t/g

:%s/\n\s*&\s*-/ -\r\t/g
:%s/\n\s*-/ -\r\t/g

"Change ** to ^
:%s/\*\*/\^/g

"Remove trailing then's from if loops
:%s/then$//g

"Replace endif, enddo etc
:%s/enddo/end/gi
:%s/endif/end/gi
"Replace .gt. .lt. etc
:%s/\.gt\./>/gi
:%s/\.lt\./</gi
:%s/\.eq\./==/gi
:%s/\.le\./<=/gi
:%s/\.ge\./>=/gi
:%s/\.true\./true/gi
:%s/\.false\./false/gi
"Replace do loops with for
:%s/do\s*[0-9a-zA-Z]*\s*=\s*[0-9a-zA-Z]*\zs,\ze/:/g
:%s/\zsdo\ze\s*[0-9a-zA-Z]*\s*=/for/g


"Replace () with [] specific to TASOPT
:%s-\parg(\([^)]\+\))-parg[\1]-g
:%s-\pari(\([^)]\+\))-pari[\1]-g
:%s-\para(\([^)]\+\))-para[\1]-g
:%s-\parm(\([^)]\+\))-parm[\1]-g
:%s-\pare(\([^)]\+\))-pare[\1]-g
