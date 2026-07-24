/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the test suite of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or (at your option) the GNU General
** Public license version 3 or any later version approved by the KDE Free
** Qt Foundation. The licenses are as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-2.0.html and
** https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

import QtQuick 2.0
import QtTest 1.1

/*!
    \qmltype SignalSpy
    \inqmlmodule QtTest
    \brief Enables introspection of signal emission.
    \since 4.8
    \ingroup qtquicktest

    In the following example, a SignalSpy is installed to watch the
    "clicked" signal on a user-defined Button type.  When the signal
    is emitted, the \l count property on the spy will be increased.

    \code
    Button {
        id: button
        SignalSpy {
            id: spy
            target: button
            signalName: "clicked"
        }
        TestCase {
            name: "ButtonClick"
            function test_click() {
                compare(spy.count, 0)
                button.clicked();
                compare(spy.count, 1)
            }
        }
    }
    \endcode

    The above style of test is suitable for signals that are emitted
    synchronously.  For asynchronous signals, the wait() method can be
    used to block the test until the signal occurs (or a timeout expires).

    \sa {QtTest::TestCase}{TestCase}, {Qt Quick Test}
*/

Item {
    id: spy
    visible: false

    TestUtil {
        id: util
    }
    // Public API.
    /*!
        \qmlproperty object SignalSpy::target

        This property defines the target object that will be used to
        listen for emissions of the \l signalName signal.

        \sa signalName, count
    */
    property var target: null
    /*!
        \qmlproperty string SignalSpy::signalName

        This property defines the name of the signal on \l target to
        listen for.

        \sa target, count
    */
    property string signalName: ""
    /*!
        \qmlproperty int SignalSpy::count

        This property defines the number of times that \l signalName has
        been emitted from \l target since the last call to clear().

        \sa target, signalName, clear()
        \readonly
    */
    readonly property alias count: spy.qtest_count
    /*!
        \qmlproperty bool SignalSpy::valid

        This property defines the current signal connection status. It will be true when the \l signalName of the \l target is connected successfully, otherwise it will be false.

        \sa count, target, signalName, clear()
        \readonly
    */
    readonly property alias valid:spy.qtest_valid
    /*!
        \qmlproperty list SignalSpy::signalArguments

        This property holds a list of emitted signal arguments. Each emission of the signal will append one item to the list, containing the arguments of the signal.
        When connecting to a new \l target or new \l signalName or calling the \l clear() method, the \l signalArguments will be reset to empty.

        \sa signalName, clear()
        \readonly
    */
    readonly property alias signalArguments:spy.qtest_signalArguments

    /*!
        \qmlmethod SignalSpy::clear()

        Clears \l count to 0, resets \l valid to false and clears the \l signalArguments to empty.

        \sa count, wait()
    */
    function clear() {
        qtest_count = 0
        qtest_expectedCount = 0
        qtest_signalArguments = []
    }

    /*!
        \qmlmethod SignalSpy::wait(timeout = 5000)

        Waits for the signal \l signalName on \l target to be emitted,
        for up to \a timeout milliseconds.  The test case will fail if
        the signal is not emitted.

        \code
        SignalSpy {
            id: spy
            target: button
            signalName: "clicked"
        }

        function test_async_click() {
            ...
            // do something that will cause clicked() to be emitted
            ...
            spy.wait()
            compare(spy.count, 1)
        }
        \endcode

        There are two possible scenarios: the signal has already been
        emitted when wait() is called, or the signal has not yet been
        emitted.  The wait() function handles the first scenario by immediately
        returning if the signal has already occurred.

        The clear() method can be used to discard information about signals
        that have already occurred to synchronize wait() with future signal
        emissions.

        \sa clear(), TestCase::tryCompare()
    */
    function wait(timeout) {
        if (timeout === undefined)
            timeout = 5000
        var expected = ++qtest_expectedCount
        var i = 0
        while (i < timeout && qtest_count < expected) {
            qtest_results.wait(50)
            i += 50
        }
        var success = (qtest_count >= expected)
        if (!qtest_results.verify(success, "wait for signal " + signalName, util.callerFile(), util.callerLine()))
            throw new Error("QtQuickTest::fail")
    }

    // Internal implementation detail follows.

    TestResult { id: qtest_results }

    onTargetChanged: {
        qtest_update()
    }
    onSignalNameChanged: {
        qtest_update()
    }

    /*! \internal */
    property var qtest_prevTarget: null
    /*! \internal */
    property string qtest_prevSignalName: ""
    /*! \internal */
    property int qtest_expectedCount: 0
    /*! \internal */
    property var qtest_signalArguments:[]
    /*! \internal */
    property int qtest_count: 0
    /*! \internal */
    property bool qtest_valid:false
    /*! \internal */

    /*! \internal */
    function qtest_update() {
        if (qtest_prevTarget != null) {
            var prevHandlerName = qtest_signalHandlerName(qtest_prevSignalName)
            var prevFunc = qtest_prevTarget[prevHandlerName]
            if (prevFunc)
                prevFunc.disconnect(spy.qtest_activated)
            qtest_prevTarget = null
            qtest_prevSignalName = ""
        }
        if (target != null && signalName != "") {
            // Look for the signal name in the object
            var func = target[signalName]
            if (typeof func !== "function") {
                // If it is not a function, try looking for signal handler
                // i.e. (onSignal) this is needed for cases where there is a property
                // and a signal with the same name, e.g. Mousearea.pressed
                func = target[qtest_signalHandlerName(signalName)]
            }
            if (func === undefined) {
                spy.qtest_valid = false
                console.log("Signal '" + signalName + "' not found")
            } else {
                qtest_prevTarget = target
                qtest_prevSignalName = signalName
                func.connect(spy.qtest_activated)
                spy.qtest_valid = true
                spy.qtest_signalArguments = []
            }
        } else {
            spy.qtest_valid = false
        }
    }

    /*! \internal */
    function qtest_activated() {
        ++qtest_count
        spy.qtest_signalArguments[spy.qtest_signalArguments.length] = arguments
    }

    /*! \internal */
    function qtest_signalHandlerName(sn) {
        if (sn.substr(0, 2) === "on" && sn[2] === sn[2].toUpperCase())
            return sn
        return "on" + sn.substr(0, 1).toUpperCase() + sn.substr(1)
    }
}
