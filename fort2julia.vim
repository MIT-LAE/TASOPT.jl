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
:%s/\n\s*&\s*\*/ \*\r\t/g
:%s/\n\s*&\s*\// \/\r\t/g
:%s/\n\s*&\s*-/ -\r\t/g

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
:%s/\.leq\./<=/gi
:%s/\.geq\./>=/gi
"Replace do loops with for
:%s/do\s*[0-9a-zA-Z]*\s*=\s*[0-9a-zA-Z]*\zs,\ze/:/g
:%s/\zsdo\ze\s*[0-9a-zA-Z]*\s*=/for/g
