MODULE MEX_INTERFACES
        INTERFACE
        subroutine mexErrMsgTxt(errormsg) bind(C, name = 'mexErrMsgTxt')
            use iso_c_binding, only: c_char
            character(c_char) :: errormsg(*)
        end subroutine mexErrMsgTxt

        subroutine mexWarnMsgTxt(warningmsg) bind(C, name = 'mexWarnMsgTxt')
            use iso_c_binding, only: c_char
            character(c_char) :: warningmsg(*)
        end subroutine mexWarnMsgTxt
    END INTERFACE
END MODULE MEX_INTERFACES
