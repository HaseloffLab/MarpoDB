function showTooltip(elem, msg) {
    elem.setAttribute('class', 'button right tooltipped');
    elem.setAttribute('aria-label', msg);
}

function fallbackMessage(action) {
    var actionMsg = '';
    var actionKey = (action == 'cut' ? 'X' : 'C');
    if (/iPhone|iPad/i.test(navigator.userAgent)) {
        actionMsg = 'No support :(';
    } else if (/Mac/i.test(navigator.userAgent)) {
        actionMsg = 'Press ⌘-' + actionKey + ' to ' + action;
    } else {
        actionMsg = 'Press Ctrl-' + actionKey + ' to ' + action;
    }
    return actionMsg;
}