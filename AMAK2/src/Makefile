ifdef ComSpec
    RM=del /F /Q
    COPY=copy
else
    RM=rm -rf
    COPY=cp
endif

all: amak ms main

ms: ms.tpl 
	@admb -f ms.tpl 

main: main.tpl 
	@admb -f main.tpl 

amak: amak.tpl 
	@admb  amak.tpl 

clean:
	@$(RM) *.cpp
	@$(RM) *.htp
	@$(RM) *.obj
